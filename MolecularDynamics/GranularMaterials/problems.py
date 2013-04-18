#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
# import matplotlib as pl
import math
# from scipy.spatial import KDTree
##import unittest
##import itertools
##from pprint import pprint, pformat
# from collections import defaultdict#, namedtuple
#from matplotlib.animation import FuncAnimation  # v1.1+
#import timeit
#import cPickle as pickle
import os

from patku_sim import Container, moldyn#, graphical, VerletIntegrator,
from patku_sim.libs import Struct, distance, make_dirs_if_necessary


class SimStats(object):
    def __init__(self, info_for_naming, containers, sim_wide_params,
            **attributes):
        self.info_for_naming = info_for_naming
        self.sim_wide_params = sim_wide_params
        self.__dict__.update(attributes)
        # Set up multi-dim data columns
        c0 = containers[0]
        y_hole = sim_wide_params.y_funnel_bottom
        y_bottom_pbc = c0.Ly0
##        flux_height = distance(y_bottom_pbc, sim_wide_params.y_funnel_bottom)

        self.cols_1d = ['times', 'grains_below_aperture']#, 'flux']
        # Some of the following columns are n-dimensional (e.g. accels would have a column for each anchor)
        self.cols_nd = ['anchor_accels']  # n-dim
        self.attr_names = self.cols_1d + self.cols_nd
        for attr in self.attr_names:
            setattr(self, attr, [])
        for c in containers[1:]:  # first container usually has no useful information (and is missing attributes)
            self.times.append(c.time)
            self.anchor_accels.append(c.anchor_accels)
            grains_below_aperture = sum(c.positions < y_hole)[1]  # only in y
            self.grains_below_aperture.append(grains_below_aperture)
##            self.flux.append(grains_below_aperture / flux_height)
        for attr in self.attr_names:
            # convert to arrays
            setattr(self, attr, np.array(getattr(self, attr)))
        # Flux stuff
        self.cumulative_grains_below_aperture = np.cumsum(self.grains_below_aperture)


    def write_csvs(self, take_every=1):
        pre = 'dumps/' + self.info_for_naming + '/'  # path prefix
        make_dirs_if_necessary(pre)
        for func, fname in zip(
                [self.anchor_csv_iter, self.flux_csv_iter],
                ['anchor_accels.csv', 'grain_count.csv']):
            with open(pre + fname, 'w') as f:
                for line in func(take_every):
                    f.write(line + '\n')

    def flux_csv_iter(self, take_every=1):
        print 'fluxish hist', np.histogram(self.cumulative_grains_below_aperture)
        yield ','.join(['time', 'cumulative_grains_below_aperture'])
        for time_ix, time in enumerate(self.times):
            if time_ix % take_every != 0:
                continue
            yield ','.join(str(a) for a in [time, self.cumulative_grains_below_aperture[time_ix]])

    def anchor_csv_iter(self, take_every=1):
##        'time', 'ix', 'value'
##        0, 0, 5
##        0, 1, 2
##        ...
##        9, 1, 2
        yield ','.join(['time', 'ix', 'value'])
        for time_ix, time in enumerate(self.times):
            if time_ix % take_every != 0:
                continue
            for ix, value in enumerate(self.anchor_accels[time_ix]):
                yield ','.join(str(a) for a in [time, ix, value])

    def __str__(self):
        return self.info_for_naming

    def __repr__(self):
        keep = self.__dict__.viewkeys() - {'times', 'containers', 'info_for_naming'}
        return '{} {}'.format(self.info_for_naming,
                {k:v for k,v in self.__dict__.iteritems() if k in keep})


class GravityForce(moldyn.Force):
    name = 'Gravity'
    def __init__(self, g=3):
        self.g = g

    def __call__(self, container, sim_wide_params):
        # Only affect y acceleration
        axs = np.zeros(container.num_particles)
        ays = np.repeat(-self.g, container.num_particles)
        if len(container.dims_iter()) == 3:
            azs = axs
            all_as = np.array([axs, ays, azs]).transpose()
        else:
            all_as = np.array([axs, ays]).transpose()
        all_as[sim_wide_params.anchor_ixs] = 0
        return all_as


class LJRepulsiveAndDampingForce(moldyn.Force):
    name = 'LJ Repulsive and Damping'
    attrs_to_copy = ['pe', 'anchor_ixs']

    def __init__(self, cutoff_dist, radius=1.0, epsilon=1.0, viscous_damping_magnitude=0):
##        super(type(self), self).__init__()
        self.cutoff_dist = cutoff_dist
        self.radius = radius
        self.epsilon = epsilon
        self.damping_magnitude = viscous_damping_magnitude

    def __call__(self, container, sim_wide_params):
        distance_matrices = container.get_distance_matrices()
        # LJ
        accelerations, pe = moldyn.lennardJonesForce(
                distance_matrices,
                radius=self.radius,
                epsilon=self.epsilon,
                cutoff_dist=self.cutoff_dist,
                anchor_ixs = sim_wide_params.anchor_ixs)
        self.pe = pe  # should LJ func calc this differently since we are doing a cutoff?
        # Damping
        dx, dy, dr = distance_matrices  # KLUDGE
        ixs_not_touching = dr > self.cutoff_dist
        self.ixs_not_touching = ixs_not_touching
        dvx, dvy, dvr = moldyn.get_distance_matrices(container.velocities)
##        if np.any(container.velocities > 1):
##            pass
##        if container.time > 1:
##            pass
        dr[dr == 0] = 0.1  # prevent / 0 errors. Will this cause problems?
        accel = self.damping_magnitude * (dvx*dx + dvy*dy)/dr
        accel[ixs_not_touching] = 0
        accelerations[:, 0] += np.sum(np.triu(accel * dx/dr), axis=0)
        accelerations[:, 1] += np.sum(np.triu(accel * dy/dr), axis=0)
        return accelerations


def hourglass(funnel_width, hole_width, d_angle_from_horizon_to_wall, dist_between_anchors, diam, **kwargs_ignore):
    """ Make a container with "anchor" particles set up as walls of an hourglass
    :param kwargs_ignore: ignored extra args so you can use an exploded dictionary (**myargdict) for arg input
    :param diam: diameter of each grain and anchor
    :returns: container, y_funnel_bottom
    """
    r_angle = math.radians(d_angle_from_horizon_to_wall)
    tan_over_2 = math.tan(r_angle)/2.0
    height_funnel = funnel_width * tan_over_2
    hole_width += diam  # trick to allow space so bottom anchors don't overlap, e.g. so a zero-width hole actually places anchors with their borders touching instead of with their centers touching
    height_hole = hole_width * tan_over_2
    xdist = dist_between_anchors * math.cos(r_angle)
    ydist = dist_between_anchors * math.sin(r_angle)
    # Centers of wall particles
    def get_anchor_centers(xstart, xend, xdist, left_side=True):
        # np.arange doesn't let you go from a greater to a lesser value, so invert twice to get that behavior
        if xstart > xend:
            dx = (xstart - xend)
##            xs = -np.arange(-xstart, -xend, xdist)
        else:
            dx = (xend - xstart)
##            xs = np.arange(xstart, xend, xdist)
        num_anchors = math.ceil(dx / xdist)
        xs = np.linspace(xstart, xend, num_anchors)
        return xs
    # Right wall equation: goes from ( (wf + wh)/2, hh - hf ) to (wf, 0)
    #m = (height_funnel - height_hole) / (funnel_width - (funnel_width + hole_width)/2.0)
    extra_funnel_extension = 0  # messes up hole alignment
    y_funnel_bottom = height_hole - height_funnel  # y is negative
    xrs = get_anchor_centers((funnel_width + hole_width)/2.0, funnel_width + extra_funnel_extension, xdist)
    yrs = get_anchor_centers(y_funnel_bottom, extra_funnel_extension, ydist, False)
    # Left wall goes from (0, 0) to ( (wf - wh)/2, hh - hf)
    xls = get_anchor_centers(-extra_funnel_extension, (funnel_width - hole_width)/2.0, xdist)
    yls = get_anchor_centers(extra_funnel_extension, y_funnel_bottom, ydist, False)

    # Place anchors
    c = Container()
    def add_anchors(xs, ys):
        for x, y in zip(xs, ys):
            c.add_particle([x, y])
    add_anchors(xls, yls)
    add_anchors(xrs, yrs)
    return c, y_funnel_bottom


def add_sand(container, width, height, dx, dy):
    moldyn.add_triangle_lattice(container, (0, width), (1, height+1), dx, dy)
    return container


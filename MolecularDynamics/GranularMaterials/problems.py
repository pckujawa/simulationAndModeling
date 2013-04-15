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

from patku_sim import Container, moldyn#, graphical, VerletIntegrator,
from patku_sim.libs import Struct


class GravityForce(moldyn.Force):
    name = 'Gravity'
    def __init__(self, g=3):
        self.g = g

    def __call__(self, container):
        # Only affect y acceleration
        axs = np.zeros(container.num_particles)
        ays = np.repeat(-self.g, container.num_particles)
        if len(container.dims_iter()) == 3:
            azs = axs
            return np.array([axs, ays, azs]).transpose()
        return np.array([axs, ays]).transpose()


class LJRepulsiveAndDampingForce(moldyn.Force):
    name = 'LJ Repulsive and Damping'
    attrs_to_copy = ['pe', 'anchor_ixs']

    def __init__(self, cutoff_dist, radius=1.0, epsilon=1.0, viscous_damping_magnitude=0, anchor_ixs=None):
##        super(type(self), self).__init__()
        self.cutoff_dist = cutoff_dist
        self.radius = radius
        self.epsilon = epsilon
        self.damping_magnitude = viscous_damping_magnitude
        self.anchor_ixs = anchor_ixs or []

    def __call__(self, container):
        distance_matrices = container.get_distance_matrices()
        # LJ
        accelerations, pe = moldyn.lennardJonesForce(
            distance_matrices, radius=self.radius, epsilon=self.epsilon, cutoff_dist=self.cutoff_dist, anchor_ixs = self.anchor_ixs)
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


def hourglass(funnel_width, hole_width, d_angle_from_horizon_to_wall, dist_between_anchors, **kwargs_ignore):
    """ Make a container with "anchor" particles set up as walls of an hourglass
    :param: kwargs_ignore ignored extra args so you can use an exploded dictionary (**myargdict) for arg input
    :returns: container, y_funnel_bottom
    """
    r_angle = math.radians(d_angle_from_horizon_to_wall)
    tan_over_2 = math.tan(r_angle)/2.0
    height_funnel = funnel_width * tan_over_2
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

    # Boundary so particles keep recirculating
    #  x is just a fudge value > funnel width so it doesn't affect anything
    #  y is some arbitrary line below the bottom of the funnel
##    c.bounds = (funnel_width+1, y_funnel_bottom - 5*hole_width)
    #TODO check that bounds code will allow negative y
    return c, y_funnel_bottom

def add_sand(container, width, height, dx, dy):
    moldyn.add_triangle_lattice(container, (0, width), (1, height+1), dx, dy)
    return container


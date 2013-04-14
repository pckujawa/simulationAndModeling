#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib as pl
import math
from scipy.spatial import KDTree
##import unittest
##import itertools
##from pprint import pprint, pformat
from collections import defaultdict#, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+
import timeit
import cPickle as pickle

from patku_sim import Container, VerletIntegrator, moldyn, graphical
from patku_sim.libs import Struct
import problems

nice_image_codecs = ['libx264']  # keep everything sharp and aren't slow
anim_save_kwargs = {'fps': 30, 'codec': nice_image_codecs[0]}
also_run_backwards = False
show_animation = True
save_animation = False
particle_radius = 0.5*2**(1.0/6)
diam = 2*particle_radius
frame_show_modulus = 10  # only show every nth frame
figsize = (4, 10)  # for animation only
num_frames_to_bootstrap = 1

dt = 1e-2
p = Struct(name='Sand Sim parameters',
    gravity_magnitude = 3,
    funnel_width = 20*diam,
    hole_width = 4*diam,
    d_angle_from_horizon_to_wall = 45,  # d for degrees
    dist_between_anchors = diam,
    num_grains = 50)  # may use dx, dy for tri lattice instead of N

xlim = (0, p.funnel_width)
ylim = tuple(np.array([0, p.num_grains]) - p.num_grains/2.0)
lstats = []  # list of stats



class SimStats(object):
    def __init__(self, info_for_naming, containers, **attributes):
        self.info_for_naming = info_for_naming
        self.__dict__.update(attributes)
        times = []
        for c in containers[1:]:  # first container usually has no useful information (and is missing attributes)
            times.append(c.time)
        self.times = np.array(times)

    def __str__(self):
        return self.info_for_naming

    def __repr__(self):
        keep = self.__dict__.viewkeys() - {'times', 'containers', 'info_for_naming'}
        return '{} {}'.format(self.info_for_naming,
                {k:v for k,v in self.__dict__.iteritems() if k in keep})


class GravityForce(object):
    def __init__(self, g=3):
        self.g = g

    def __call__(self, container):
        # Only affect y acceleration
        axs = np.zeros(container.num_particles)
        ays = np.repeat(g, container.num_particles)
        azs = axs
        return axs, ays, azs


class RunFunc():
    def __init__(self):
        pass
    def __call__(self):
        return False


def run(
    normal_force = 0,
    neighbor_facilitator = None
    ):
    global containers, lstats
    lj_method = 'dm' if neighbor_facilitator is None else 'nl'  # LJ method used either distance matrix or neighbor list
    info_for_naming = '{} g={:.0f} dt={:.2f}'.format('numParticles', normal_force, dt)

    print 'running with', info_for_naming

    init_container, y_funnel_bottom = problems.hourglass(**p.__dict__)
    problems.add_sand(init_container, p.funnel_width, diam, diam)
    containers = [init_container]
    integrator = VerletIntegrator()

    run_func = None# RunFunc(last_particle_position)
    graphical.animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation, show_animation, run_func,
    anim_save_kwargs = anim_save_kwargs)

    stats = SimStats(info_for_naming, containers,
        gravity_magnitude = normal_force,
        dt = dt)

    lstats.append(stats)
    # graphical.plot_sand(stats.times, stats.???, normal_force=stats.gravity_magnitude, info_for_naming=info_for_naming, show=False)


for use_neighbors in [False]:
    neighbor_facilitator = NeighborFacilitator(2.5, 10) if use_neighbors else None
    run_timed = lambda: run(p.gravity_magnitude, neighbor_facilitator)
    num_times = 1
    timer = timeit.Timer(run_timed)
    time_taken = timer.timeit(num_times)
    used_ns = 'with' if use_neighbors else 'without'
    print '\tfinished {} timed run {} neighbors in {}s'.format(num_times, used_ns, time_taken)

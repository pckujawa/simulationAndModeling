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
save_animation = True
particle_radius = 0.5*2**(1.0/6)
diam = 2*particle_radius
frame_show_modulus = 10  # only show every nth frame
figsize = (6, 10)  # for animation only
num_frames_to_bootstrap = 1

dt = 1e-2
p = Struct(name='Sand Sim parameters',
    dt = dt,  # for easier naming
    gravity_magnitude = 3,
    viscous_damping_magnitude = -10,
    funnel_width = 10*diam,
    hole_width = diam,#3*diam,  # can't be zero or particles overlap
    d_angle_from_horizon_to_wall = 60,  # d for degrees
    dist_between_anchors = diam+0.01,
    num_grains = 100)  # may use dx, dy for tri lattice instead of N

def params_tostring(self):
    return '{num_grains} grains g={gravity_magnitude:.1f} damp={viscous_damping_magnitude} hw={hole_width:.2f} angle={d_angle_from_horizon_to_wall} dt={dt:.2f}'.format(**self.__dict__)

xlim = (-diam, p.funnel_width + diam)
ylim = np.array([-1.1, 0.9]) * p.funnel_width/2.0 * math.tan(math.radians(p.d_angle_from_horizon_to_wall))
lstats = []  # list of stats


class RunFunc():
    def __init__(self, y_funnel_bottom):
        self.y_funnel_bottom = y_funnel_bottom

    def __call__(self, container):
        return False
        count = np.sum(container.positions < self.y_funnel_bottom)
        print count, 'grains thru the hole'
        return count < p.num_grains



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




def run(
    neighbor_facilitator = None
    ):
    global containers, lstats, ixs_unmoving
    lj_method = 'dm' if neighbor_facilitator is None else 'nl'  # LJ method used either distance matrix or neighbor list
    info_for_naming = params_tostring(p)

    print 'running with', info_for_naming

    init_container, y_funnel_bottom = problems.hourglass(**p.__dict__)
    ixs_unmoving = range(init_container.num_particles)
    forces = [
        problems.GravityForce(g = p.gravity_magnitude),
        problems.LJRepulsiveAndDampingForce(cutoff_dist = diam,
                viscous_damping_magnitude = p.viscous_damping_magnitude,
                anchor_ixs = ixs_unmoving)
    ]

    packing_factor = 1.0  # guess as to how many particles fit into a unit area
    sand_height = p.num_grains / (p.funnel_width * packing_factor)
    problems.add_sand(init_container, p.funnel_width, sand_height, diam, diam)
    containers = [init_container]

    integrator = VerletIntegrator()
    integrator.forces = forces
    integrator.ixs_unmoving = ixs_unmoving

    run_func = RunFunc(y_funnel_bottom)

    sim_wide_params = Struct(
            anchor_ixs = ixs_unmoving)

    graphical.animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation, show_animation, run_func,
        anim_save_kwargs = anim_save_kwargs,
        sim_wide_params = sim_wide_params)

    stats = SimStats(info_for_naming, containers,
        gravity_magnitude = p.gravity_magnitude,
        dt = dt)

    lstats.append(stats)
    # graphical.plot_sand(stats.times, stats.???, gravity_magnitude=stats.gravity_magnitude, info_for_naming=info_for_naming, show=False)


for use_neighbors in [False]:
    neighbor_facilitator = NeighborFacilitator(2.5, 10) if use_neighbors else None
    run_timed = lambda: run(neighbor_facilitator)
    num_times = 1
    timer = timeit.Timer(run_timed)
    try:
        time_taken = timer.timeit(num_times)
        used_ns = 'with' if use_neighbors else 'without'
        print '\tfinished {} timed run {} neighbors in {}s'.format(num_times, used_ns, time_taken)
    except AttributeError as e:
        print e  # TkInter error usually; can be ignored
        print "Can't save when in a Tk python session, but we get this error"

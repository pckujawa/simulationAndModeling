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
show_animation = False
save_animation = False  # can't show and save, for some reason (TkError)
particle_radius = 0.5*2**(1.0/6)
diam = 2*particle_radius
frame_show_modulus = 2  # only show every nth frame
figsize = (6, 10)  # for animation only
num_frames_to_bootstrap = 1

dt = 1e-2
p = Struct(name='Sand Sim parameters',
    # repeats are stored here for easier access inside various functions
    dt = dt,
    diam = diam,
    data_dump_modulus = 10,  # only save every nth container's value
    gravity_magnitude = 3,
    viscous_damping_magnitude = -10,
    funnel_width = 15*diam,
    hole_width = 0*diam,
    d_angle_from_horizon_to_wall = 60,  # d for degrees
    dist_between_anchors = diam,
    num_grains = 100,  # may use dx, dy for tri lattice instead of N
    top_bound = 10*diam  # need to be above stacks grains or they will wrap below
)

def params_tostring(self):
    hw = self.hole_width / self.diam  # print hole width proportional to d
    return 'hw={hw:.2f}d angle={d_angle_from_horizon_to_wall} grains={num_grains} g={gravity_magnitude:.1f} damp={viscous_damping_magnitude} dt={dt:.2f}'.format(hw=hw, **self.__dict__)

xlim = (-diam, p.funnel_width + diam)
ylim = np.array([-1.1, 0.9]) * p.funnel_width/2.0 * math.tan(math.radians(p.d_angle_from_horizon_to_wall))
lstats = []  # list of stats


class RunFunc():
    def __init__(self, sim_wide_params):
        self.sim_wide_params = sim_wide_params

    def __call__(self, container):
##        return False
        return container.time < 20
        count = np.sum(container.positions < self.sim_wide_params.y_funnel_bottom)
##        print count, 'grains thru the hole'
        return count < p.num_grains and container.time < 50  # failsafe




def run(neighbor_facilitator = None):
    global containers, lstats, p
    lj_method = 'dm' if neighbor_facilitator is None else 'nl'  # LJ method used either distance matrix or neighbor list
    info_for_naming = params_tostring(p)

    print 'running with', info_for_naming

    init_container, y_funnel_bottom = problems.hourglass(**p.__dict__)
    print 'y_funnel_bottom =', y_funnel_bottom
    anchor_ixs = range(init_container.num_particles)
    forces = [
        problems.GravityForce(g = p.gravity_magnitude),
        problems.LJRepulsiveAndDampingForce(cutoff_dist = diam,
                viscous_damping_magnitude = p.viscous_damping_magnitude)
    ]

    packing_factor = 1.0  # guess as to how many particles fit into a unit area
    sand_height = p.num_grains / (p.funnel_width * packing_factor)
    problems.add_sand(init_container, p.funnel_width, sand_height, diam, diam)
    containers = [init_container]

    p.anchor_ixs = anchor_ixs
    p.y_funnel_bottom = y_funnel_bottom

    integrator = VerletIntegrator()
    integrator.sim_wide_params = p
    integrator.forces = forces

    run_func = RunFunc(p)

    try:
        graphical.animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation, show_animation, run_func,
            anim_save_kwargs = anim_save_kwargs,
            sim_wide_params = p)
    except AttributeError as e:
        print e
        print "Can't save animation after showing it when in a Tk python session, but we get this error"

    stats = problems.SimStats(info_for_naming, containers, p)
    lstats.append(stats)

    stats.write_csvs(p.data_dump_modulus)
    # graphical.plot_sand(stats.times, stats.???, gravity_magnitude=stats.gravity_magnitude, info_for_naming=info_for_naming, show=False)


for use_neighbors in [False]:
    neighbor_facilitator = NeighborFacilitator(2.5, 10) if use_neighbors else None
    run_timed = lambda: run(neighbor_facilitator)
    num_times = 1
    timer = timeit.Timer(run_timed)
    time_taken = timer.timeit(num_times)
    used_ns = 'with' if use_neighbors else 'without'
    print '\tfinished {} timed run {} neighbors in {}s'.format(num_times, used_ns, time_taken)
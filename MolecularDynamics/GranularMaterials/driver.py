#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib as pl
import math
##from scipy.spatial import KDTree
##import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict#, namedtuple
##from matplotlib.animation import FuncAnimation  # v1.1+
import timeit
import cPickle as pickle

from patku_sim import Container, VerletIntegrator, moldyn, graphical
from patku_sim.libs import Struct, distance
import problems

nice_image_codecs = ['libx264']  # keep everything sharp and aren't slow
anim_save_kwargs = {'fps': 30, 'codec': nice_image_codecs[0]}
show_animation = False
save_animation = True  # can't show and save, for some reason (TkError)
particle_radius = 0.5*2**(1.0/6)
diam = 2*particle_radius
frame_show_modulus = 5  # only show every nth frame
figheight = 6  # for animation only
num_frames_to_bootstrap = 1

dt = 1e-2
p = Struct(name='Sand Sim parameters',
    # repeats are stored here for easier access inside various functions
    dt = dt,
    diam = diam,
    dump_data = True,
    data_dump_modulus = 10,  # only save every nth container's value
    gravity_magnitude = 3,
    viscous_damping_magnitude = -10,
    funnel_width = 10*diam,  # leads to odd spacing for various angles
    hole_width = [0, 1.5*diam],
    d_angle_from_horizon_to_wall = [15, 30, 45, 60],  # d for degrees
    dist_between_anchors = diam,
    num_grains = 100  # not exact though, due to packing
)

class RunFunc():
    def __init__(self, sim_wide_params):
        self.sim_wide_params = sim_wide_params

    def __call__(self, container):
##        return False
        return container.time < 1
        count = np.sum(container.positions < self.sim_wide_params.y_funnel_bottom)
##        print count, 'grains thru the hole'
        return count < p.num_grains and container.time < 50  # failsafe


def params_tostring(self):
    hw = self.hole_width / self.diam  # print hole width proportional to d
    return 'hw={hw:.2f}d angle={d_angle_from_horizon_to_wall} grains={num_grains} g={gravity_magnitude:.1f} damp={viscous_damping_magnitude} dt={dt:.2f}'.format(hw=hw, **self.__dict__)

lstats = []  # list of stats





def run(neighbor_facilitator = None):
    global containers, lstats, p
    lj_method = 'dm' if neighbor_facilitator is None else 'nl'  # LJ method used either distance matrix or neighbor list
    info_for_naming = params_tostring(p)

    print 'running with', info_for_naming

    init_container, p.y_funnel_bottom = problems.hourglass(**p.__dict__)

    p.anchor_ixs = range(init_container.num_particles)

    packing_factor = 1.25  # educated guess as to how many particles fit into a unit area
    p.sand_height = p.num_grains * packing_factor / p.funnel_width
    problems.add_sand(init_container, p.funnel_width, p.sand_height, p.diam, p.diam)
    print 'num grains expected/actual =', p.num_grains, '/', (init_container.num_particles - len(p.anchor_ixs))

    # Boundary so particles keep recirculating
    # Need to have enough x space that grains and anchors don't init on top of each other
    init_container.bounds = [(-0.1, p.funnel_width + 0.1),
            (p.y_funnel_bottom - p.diam, p.sand_height + p.diam)]

    forces = [
        problems.GravityForce(g = p.gravity_magnitude),
        problems.LJRepulsiveAndDampingForce(cutoff_dist = p.diam,
                viscous_damping_magnitude = p.viscous_damping_magnitude)
    ]

    xlim, ylim = init_container.bounds
    w = distance(*xlim)
    h = distance(*ylim)
    a_ratio = w / h
    # figsize must be integers or animation will error saving to file!
    figheight = int(figheight)
    figsize = (int(a_ratio * figheight)+1, figheight)

    containers = [init_container]

    integrator = VerletIntegrator()
    integrator.sim_wide_params = p
    integrator.forces = forces

    run_func = RunFunc(p)

    try:
        graphical.animate_with_live_integration(
            containers, integrator, dt, xlim, ylim, figsize,
            particle_radius, frame_show_modulus,
            num_frames_to_bootstrap, info_for_naming,
            save_animation, show_animation, run_func,
            anim_save_kwargs = anim_save_kwargs,
            sim_wide_params = p)
    except AttributeError as e:
        print e
        print "Can't save animation after showing it when in a Tk python session, but we get this error"

    stats = problems.SimStats(info_for_naming, containers, p)
    lstats.append(stats)

    if p.dump_data:
        stats.write_csvs(p.data_dump_modulus)
    # graphical.plot_sand(stats.times, stats.???, gravity_magnitude=stats.gravity_magnitude, info_for_naming=info_for_naming, show=False)

    print 'all params:', p


def main():
    run_timed = lambda: run()
    num_times = 1
    timer = timeit.Timer(run_timed)
    time_taken = timer.timeit(num_times)
    print '\tfinished in {t}s'.format(t = time_taken)
##    used_ns = 'with' if use_neighbors else 'without'
##    print '\tfinished {n} timed run {with_} neighbors in {t}s'.format(
##            n = num_times, with_ = used_ns, t = time_taken)

angles = p.d_angle_from_horizon_to_wall
holes = p.hole_width
try:
    for angle in angles:
        p.d_angle_from_horizon_to_wall = angle
        try:
            for hw in holes:
                p.hole_width = hw
                main()
        except TypeError:
            main()
except TypeError:
    main()

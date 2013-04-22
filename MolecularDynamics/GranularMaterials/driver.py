#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import matplotlib.pylab as pl
import math
##from scipy.spatial import KDTree
##import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict#, namedtuple
##from matplotlib.animation import FuncAnimation  # v1.1+
import timeit
import time
import cPickle as pickle

from patku_sim import Container, VerletIntegrator, moldyn, graphical
from patku_sim.libs import Struct, distance
import problems

nice_image_codecs = ['libx264']  # keep everything sharp and aren't slow
anim_save_kwargs = {'fps': 30, 'codec': nice_image_codecs[0]}
show_animation = False
save_animation = True  # can't show and save, for some reason
particle_radius = 0.5*2**(1.0/6)
diam = 2*particle_radius
frame_show_modulus = 5  # only show every nth frame
figheight = 8  # for animation only
num_frames_to_bootstrap = 1

dt = 1e-2
p = Struct(name='Sand Sim parameters',
    # repeats are stored here for easier access inside various functions
    dt = dt,
    diam = diam,
    dump_data = not show_animation,
    use_bounds = True,
    data_dump_modulus = 20,  # only save every nth container's value
    gravity_magnitude = 2,
    viscous_damping_magnitude = -10,
    funnel_width = 11*diam,  # leads to odd spacing for various angles
    hole_width = np.array([3]) * diam,  # ===============
    d_angle_from_horizon_to_wall = [60],  # d for degrees
    dist_between_anchors = diam,
    grain_height = 20  # 20 may be best all around  # 10 is good with fw=10d hw=1.5d angle=60
##    num_grains = 10  # not exact though, due to packing
)

class RunFunc():
    def __init__(self, sim_wide_params):
        self.sim_wide_params = sim_wide_params

    def __call__(self, container):
        if show_animation:
            return False  # don't cache any frames, just show live
        p = self.sim_wide_params
        # Avoid container not having correct attributes by just always skipping a few time steps
        if container.time < dt:
            return True
        if p.use_bounds:
            return container.time < 150
        if container.cumulative_grains_below_aperture < p.num_grains:
            self.checkpoint_time = container.time
            return container.time < 120  # failsafe in case of blockage
        else:
            # at this point, checkpoint will never be updated again
            return container.time < self.checkpoint_time + 1



def params_tostring(self):
    hw = self.hole_width / self.diam  # print hole width proportional to d
    bounded = 'PBC ' if p.use_bounds else ''
    return 'hw={hw:.1f}d grain_h={grain_height} g={gravity_magnitude:.1f} damp={viscous_damping_magnitude} dt={dt:.2f}/{pbc}angle={d_angle_from_horizon_to_wall}'.format(hw=hw, pbc=bounded, **self.__dict__)

lstats = []  # list of stats



def time_now():
    return time.strftime("%a %b %d %H:%M:%S")

def run(neighbor_facilitator = None):
    global containers, lstats, p
    lj_method = 'dm' if neighbor_facilitator is None else 'nl'  # LJ method used either distance matrix or neighbor list
    info_for_naming = params_tostring(p)

    header = "Starting at {t}".format(t=time_now())
    print header
    print '-'*len(header)
    print 'running with', info_for_naming

    # Get the initial problem container and meta
    init_container, p.y_funnel_bottom, p.anchor_ixs = \
            problems.hourglass(**p.__dict__)
    p.num_grains = init_container.num_particles - len(p.anchor_ixs)


    # Boundary so particles keep recirculating
    # Need to have enough x space that grains and anchors don't init on top of each other
    posns = init_container.positions
    xl, yb = np.min(posns, axis=0)  # left, bottom
    xr, yt = np.max(posns, axis=0)  # right, top
    xlim, ylim = [(xl - p.diam, xr + p.diam),
                (yb - 2*p.diam, yt + p.diam)]
    if p.use_bounds:
        init_container.bounds = [xlim, ylim]

    forces = [
        problems.GravityForce(g = p.gravity_magnitude)
        , problems.LJRepulsiveAndDampingForce(cutoff_dist = p.diam,
                viscous_damping_magnitude = p.viscous_damping_magnitude)
        , problems.FakeStatsForce()
    ]

    w = distance(*xlim)
    h = distance(*ylim)
    a_ratio = w / h
    # figsize must be integers or animation will error saving to file!
    fh = int(figheight)
    figsize = (int(a_ratio * fh)+1, fh)

    containers = [init_container]

    integrator = VerletIntegrator()
    integrator.sim_wide_params = p
    integrator.forces = forces

    run_func = RunFunc(p)

    graphical.animate_with_live_integration(
        containers, integrator, dt, xlim, ylim, figsize,
        particle_radius, frame_show_modulus,
        num_frames_to_bootstrap, info_for_naming,
        save_animation, show_animation, run_func,
        anim_save_kwargs = anim_save_kwargs,
        sim_wide_params = p)

    stats = problems.SimStats(info_for_naming, containers, p)
    lstats.append(stats)

    if p.dump_data:
        stats.write_csvs(p.data_dump_modulus)
        stats.write_plots(p.data_dump_modulus, show=False)
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
for angle in angles:
    p.d_angle_from_horizon_to_wall = angle
    for hw in holes:
        p.hole_width = hw
        main()
##
##started = False
##try:
##    for angle in angles:
##        p.d_angle_from_horizon_to_wall = angle
##        try:
##            for hw in holes:
##                started = True
##                p.hole_width = hw
##                main()
##        except TypeError as e:# if 'iterator' in e.message:
##            print e
##            if not started:
##                started = True
##                main()
##except TypeError as e:
##    print e
##    if not started:
##        started = True
##        main()

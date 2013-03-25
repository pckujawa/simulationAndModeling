#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import pylab as pl
import math
##import unittest
##import itertools
##from pprint import pprint, pformat
from collections import defaultdict#, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+
import timeit
import cPickle as pickle

from patku_sim import Container, VerletIntegrator, moldyn
from patku_sim.libs import Struct
import graphical

nice_image_codecs = ['libx264']  # keep everything sharp and aren't slow
fast_codecs = []
bad_codecs = ['v210', 'rawvideo', '']
anim_save_kwargs = {'fps': 30, 'codec': nice_image_codecs[0]}
also_run_backwards = False
show_animation = False
save_animation = False
particle_radius = 0.5*2**(1.0/6)
a = 2*particle_radius
frame_show_modulus = 10  # only show every nth frame
dt = 1e-2
xlim, ylim = (0, 30), (-1, 5)
figsize = (10, 4)  # for animation only
sled_k = 500.0
pulling_force_k = 5.0
allow_negative_pull_force = True
damp_force_multiplier = 10.0
num_frames_to_bootstrap = 100
lstats = []  # list of stats
##unpickle_filepath = r'lstats 1_9_13_17 Fpv=0.1 W=-20..40by5.pickle'


class RunFunc():
    def __init__(self, last_particle_position):
        self.time_of_break = None
        self.last_particle_position = last_particle_position
    def __call__(self, container):
        # Run until we break static friction (move by 'a') and then some
        if self.time_of_break is not None:
            return container.time < self.time_of_break + 10
        sled_lead_posn = container.positions[-1]
        x = (sled_lead_posn - self.last_particle_position)[0]
        if x > a:
            self.time_of_break = container.time
        return True


def create(cnt_sled_particles, cnt_floor_particles=100):
    c = Container()
    # Make floor with particles spaced 'a' apart
    floor_ixs = range(cnt_floor_particles)
    for i in floor_ixs:
        c.add_particle([a*(i-1), 0])
    c.floor_particle_ixs = floor_ixs
    # Make sled with particles in equilateral triangle of side length '2a' with base of two particles placed 'a' above floor starting at a/2
    top = Struct(x = 3.0*a/2, y = a * (math.sqrt(3) + 1))
    bottom = Struct(x = a/2.0, y = a)
    for i in xrange(cnt_sled_particles):
        if i % 2 == 1:
            side = top
        else:
            side = bottom
        last_particle_position = [side.x, side.y]  # will be true when loop exits
        c.add_particle(last_particle_position)
        side.x += 2*a
    c.sled_particle_ixs = np.array(range(cnt_sled_particles)) + cnt_floor_particles
    return c, last_particle_position

#  To vary:
def run(
    normal_force = 0,
    sled_conditions = (9, 25),  # num particles in sled, floor
    pulling_force_v = 0.1
    ):
    global containers, lstats
    info_for_naming = '{} W={:.0f} Fpv={:.2f} Fpk={:.1f} dt={:.2f} k={:.1f} Fdm={:.1f}'.format(sled_conditions, normal_force, pulling_force_v, pulling_force_k, dt, sled_k, damp_force_multiplier)

    print 'running with', info_for_naming

    init_container, last_particle_position = create(*sled_conditions)
    containers = [init_container]
    integrator = VerletIntegrator()
    sled_forcer = moldyn.SledForcer(2*a, u=last_particle_position, k=sled_k)
    sled_forcer.pulling_force_v = pulling_force_v
    sled_forcer.pulling_force_k = pulling_force_k
    sled_forcer.allow_negative_pull_force = allow_negative_pull_force
    sled_forcer.damp_force_multiplier = damp_force_multiplier
    sled_forcer.normal_force = normal_force
    integrator.sled_forcer = sled_forcer

    graphical.animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation, show_animation, RunFunc(last_particle_position),
    anim_save_kwargs = anim_save_kwargs)

    stats = SimStats(info_for_naming, containers,
        sled_size = sled_conditions[0],
        W = normal_force,
        Fpv = pulling_force_v,
        Fpk = pulling_force_k,
        dt = dt,
        sled_k = sled_k,
        Fdm = damp_force_multiplier)
    lstats.append(stats)
    graphical.plot_pulling_force(stats.times, stats.pulling_forces, normal_force=stats.W, info_for_naming=info_for_naming, show=False)


class SimStats(object):
    def __init__(self, info_for_naming, containers, **attributes):
        self.info_for_naming = info_for_naming
        self.__dict__.update(attributes)
        times = []
        pulling_forces = []  # x and y
        for c in containers[1:]:  # first container usually has no useful information (and is missing attributes)
            times.append(c.time)
            pulling_forces.append(c.pull_accelerations)
        self.times = np.array(times)
        self.pulling_forces = np.array(pulling_forces)
        assert 'W' in attributes
        pfm = graphical.get_pulling_force_max_and_meta(self.times, self.pulling_forces, normal_force=self.W)
        self.Fp_max, self.Fp_max_time, ix = pfm
    def __str__(self):
        return self.info_for_naming
    def __repr__(self):
        keep = self.__dict__.viewkeys() - {'times', 'pulling_forces', 'containers', 'info_for_naming'}
        return '{} {}'.format(self.info_for_naming,
            {k:v for k,v in self.__dict__.iteritems() if k in keep})

try:
    lstats = pickle.load(open(unpickle_filepath, 'rb'))
except:
    for num_sled in [9]:# [1, 9, 13, 17]:
        for pull_v in [0.1]:# np.linspace(0.05, 0.50, 10):
            for W in [0]:# xrange(0, 41, 5):  # include 40
                run_timed = lambda: run(W, (num_sled, 25), pull_v)
                num_times = 1
                timer = timeit.Timer(run_timed)
                time_taken = timer.timeit(num_times)
                print '\tfinished {} iterations in {}s'.format(num_times, time_taken)
    pickle.dump(lstats, open('lstats.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

def get_csv_iter(lstats):
    csv_cols = ['sled_size', 'W', 'Fp_max', 'Fp_max_time', 'Fpv']
    yield ','.join(csv_cols)
    for s in lstats:
        yield ','.join([str(getattr(s, col)) for col in csv_cols])

print '\n'.join(get_csv_iter(lstats))

data_frame = np.genfromtxt(get_csv_iter(lstats), delimiter=',', names=True)
##dm = pandas.DataMatrix.fromRecords(data_frame)

if len(lstats) > 1:
    graphical.plot_friction_slope(data_frame)
graphical.plot_all_pulling_forces(lstats)

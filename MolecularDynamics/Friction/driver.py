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
##from collections import defaultdict, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+

from patku_sim import Container, VerletIntegrator, moldyn
from patku_sim.libs import Struct
import graphical

also_run_backwards = False
show_animation = False
save_animation = True
particle_radius = 2**(1.0/6)  # a.k.a. 'a'; actually, the point at which the Lenn-Jones force is stable between particles
frame_show_modulus = 5  # only show every nth frame
dt = 1e-2
sim_name = 'friction'
xlim, ylim = (0, 30), (-1, 5)
figsize = (10, 4)
sled_conditions = (9, 25)  # num particles in sled, floor
sled_k = 500.0
pulling_force_v = 0.1
pulling_force_k = 5.0
allow_negative_pull_force = True
damp_force_multiplier = 10.0
num_frames_to_bootstrap = 100
info_for_naming = '{} dt={:.2f} k={:.1f} Fpk={:.1f}(-{}) Fpv={:.1f} Fdm={:.1f}'.format(sled_conditions, dt, sled_k, pulling_force_k, allow_negative_pull_force, pulling_force_v, damp_force_multiplier)


def run_func(container):
    sled_lead_posn = container.positions[-1]
    x = (sled_lead_posn - last_particle_position)[0]
    return x < particle_radius

print 'running with', info_for_naming
last_particle_position = None
def create(cnt_sled_particles, cnt_floor_particles=100):
    """
    """
    global particle_radius, last_particle_position
    c = Container()
    a = particle_radius
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
    return c

init_container = create(*sled_conditions)
containers = [init_container]
integrator = VerletIntegrator()
sled_forcer = moldyn.SledForcer(2*particle_radius, u=last_particle_position, k=sled_k)
sled_forcer.pulling_force_v = pulling_force_v
sled_forcer.pulling_force_k = pulling_force_k
sled_forcer.allow_negative_pull_force = allow_negative_pull_force
sled_forcer.damp_force_multiplier = damp_force_multiplier
integrator.sled_forcer = sled_forcer

graphical.animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation, show_animation, run_func)

times = []
pulling_forces = []  # x and y
for c in containers[1:]:  # first container usually has no useful information (and is missing attributes)
    times.append(c.time)
    pulling_forces.append(c.pull_accelerations)

graphical.plot_pulling_force(times, pulling_forces, info_for_naming)


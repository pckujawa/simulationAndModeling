#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
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
from graphical import animate_with_live_integration, plot_potential_energy

also_run_backwards = False
show_animation = True
save_animation = True
particle_radius = 2**(1.0/6)  # a.k.a. 'a'; actually, the point at which the Lenn-Jones force is stable between particles
frame_show_modulus = 5  # only show every nth frame
dt = 1e-2
sim_name = 'friction'
xlim, ylim = (0, 30), (-1, 5) # (0, 50), (-5, 5*particle_radius)
figsize = (10, 4)
sled_conditions = (13, 25)  # sled, floor
spring_const = 5e0
pulling_force_multiplier = 5e0
allow_negative_pull_force = True
damp_force_multiplier = 1e2
num_frames_to_bootstrap = 100
info_for_naming = '{} dt={:.3f} k={:.1f} Fpm={:.1f}(-{}) Fdm={:.1f}'.format(sled_conditions, dt, spring_const, pulling_force_multiplier, allow_negative_pull_force, damp_force_multiplier)
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
sled_forcer = moldyn.SledForcer(2*particle_radius, u=last_particle_position, k=spring_const)
sled_forcer.pulling_force_multiplier = pulling_force_multiplier
sled_forcer.allow_negative_pull_force = allow_negative_pull_force
sled_forcer.damp_force_multiplier = damp_force_multiplier
integrator.sled_forcer = sled_forcer

animate_with_live_integration(containers, integrator, dt, xlim, ylim, figsize, particle_radius, frame_show_modulus, num_frames_to_bootstrap, info_for_naming, save_animation)


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

also_run_backwards = False
show_animation = True
save_animation = False
particle_radius = 1 # 2**(1.0/6)  # a.k.a. 'a'; actually, the point at which the Lenn-Jones force is stable between particles
frame_show_modulus = 10  # only show every nth frame
dt = 1e-2
sim_name = 'friction'
xlim, ylim = (0, 20), (-1, 5) # (0, 50), (-5, 5*particle_radius)
figsize = (10, 4)
sled_conditions = (5, 10)
spring_const = 10
num_frames_to_bootstrap = 100

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
        c.add_particle([a*i-1, 0])
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
    # Connect sled particles like a trianglular truss

    return c

init_container = create(*sled_conditions)
containers = [init_container]
integrator = VerletIntegrator()
integrator.sled_forcer = moldyn.SledForcer(2*particle_radius, u=last_particle_position, k=spring_const)

# Animate orbit
# Code courtesy of George Lesica
fig = pl.figure(figsize=figsize)
ax = pl.axes(xlim=xlim, ylim=ylim)  # necessary for animation to know correct bounds
ax.set_aspect('equal')
##ax.set_xlim((0, init_container.bounds[0]))
##ax.set_ylim((0, init_container.bounds[1]))
pl.title('Molec Dyn Friction Simulation', fontsize=16)
pl.xlabel('X Position')
pl.ylabel('Y Position')


## (Kevin) Pre initializing is necessary I guess
posns = init_container.positions
circles = []
for i,posn in enumerate(posns):
    e = moldyn.get_nice_circle(posn[0], posn[1], 0.5*particle_radius)
    circles.append(ax.add_patch(e))

num_forward_frames = 0
# Bootstrap some frames
while num_forward_frames < num_frames_to_bootstrap:
    num_forward_frames += 1
    next_container = integrator.step(containers[-1], dt)
    containers.append(next_container)

# init fn seems to prevent 'ghosting' of first-plotted data in others' code
def init():
    """initialize animation"""
    for c in circles:
        c.center = (-1, -1)  # hide (hopefully) off-screen
    return circles

def next_frame(ix_frame):
    global num_forward_frames
    print 'frame', ix_frame
    # Do only the integration necessary to get to the requested frame
    while num_forward_frames < ix_frame:
        num_forward_frames += 1
        for _ in xrange(1 + frame_show_modulus):  # always run at least once
            next_container = integrator.step(containers[-1], dt)
            containers.append(next_container)
    posns = containers[ix_frame].positions
    facecolor = 'green'
    # TODO paint floor black, sled green, etc
    for i,circle in zip(xrange(init_container.num_particles), circles):
        circle.center = (posns[i][0], posns[i][1])  # x and y
        circle.set_facecolor(facecolor)
    return circles

anim = FuncAnimation(fig, next_frame, interval=dt, blit=True, init_func=init)
if show_animation:
    pl.show()

end_container = containers[-1]

# Now run... backwards!  TODO
if also_run_backwards:
    for i in xrange(num_forward_frames):
        next_container = integrator.step(containers[-1], -dt)
        containers.append(next_container)


num_total_frames = num_forward_frames
if also_run_backwards:
    num_total_frames += num_forward_frames


if save_animation:
    # Seems like we need to re-run to get the full movie
    anim = FuncAnimation(fig, next_frame, interval=dt, blit=True, frames=num_total_frames)
    try:
        anim.save('mol_dyn_friction_{}_animation.avi'.format(num_forward_frames), fps=30)
    except:  # Tk error
        pass


pl.clf()
# Plot PE --------------------------
times = pl.frange(dt, num_total_frames*dt, dt)  # skip zeroth time because we have no value for it
pes = [c.potential_energy for c in containers[1:]]  # skip first container because it has no PE
plotted_pes = np.array(pes[:len(times)])
plotted_pes /= init_container.num_particles  # PE per particle
pl.plot(times, plotted_pes, 'o-', color='black', markersize=1, linewidth=0.1)
pl.ylabel('Potential energy per particle')
pl.xlabel('Time')
pl.title('PE/particle for {} frames'.format(num_forward_frames))
pl.savefig('mol_dyn_{}_pe.png'.format(num_forward_frames))
pl.show()

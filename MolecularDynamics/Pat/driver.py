#-------------------------------------------------------------------------------
# Name:        molecular dynamics driver
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
from problems import get_container_for

also_run_backwards = False
show_animation = True
save_animation = False
particle_radius = 2**(1.0/6)
num_forward_frames = 5000
frame_show_modulus = 10  # only show every nth frame
dt = 1e-2
sim_name = ['symmetry', 'problem_3', 'problem_1', 'line'][0]
extra_params = {'Lx': 10,
    'random_velocity': 0, 'random_particle_ix': 0,
    'lattice': ['triangle', 'square'][0],
    'symmetry_number': 'six'}


def squeeze(container, squeeze_factor, t):
    c = container
    Lx, Ly = c.bounds
    # c is the container object.
    if t > 3. and Lx > 8.0 * 2.0**(1.0/6):
        # Squeeze the box!
        Lx *= squeeze_factor
        Ly *= squeeze_factor
##        c.x  *=  squeeze_factor
##        c.y  *=  squeeze_factor
##where the squeeze isn't applied until the atoms settle down a little bit, and doesn't continue past the solid packing size. SQUEEZE_FACTORS > .995 work well. Only apply the squeeze about every 20 time steps
    return container

init_container, special_particles = get_container_for(sim_name, **extra_params)
print 'special_particles:', special_particles

containers = [init_container]
integrator = VerletIntegrator()
for i in xrange(num_forward_frames):
    next_container = integrator.step(containers[-1], dt)
    containers.append(next_container)
end_container = containers[-1]

# Now run... backwards!
if also_run_backwards:
    for i in xrange(num_forward_frames):
        next_container = integrator.step(containers[-1], -dt)
        containers.append(next_container)

# Animate orbit
# Code courtesy of George Lesica
fig = pl.figure(figsize=(4, 4))
xlim, ylim = init_container.bounds
ax = pl.axes(xlim=(0, xlim), ylim=(0, ylim))  # necessary because initial plot is too zoomed in
ax.set_aspect('equal')
ax.set_xlim((0, init_container.bounds[0]))
ax.set_ylim((0, init_container.bounds[1]))
pl.title('Molec Dyn Simulation', fontsize=16)
pl.xlabel('X Position')
pl.ylabel('Y Position')


## (Kevin) Pre initializing is necessary I guess
posns = init_container.positions
circles = []
for i,posn in enumerate(posns):
    e = moldyn.get_nice_circle(posn[0], posn[1], 0.5*particle_radius)
    circles.append(ax.add_patch(e))

def init():
    for c in circles:
        c.center = (-1, -1)  # hide
    return circles

def next_frame(ix_frame):
    ix_frame *= frame_show_modulus
    posns = containers[ix_frame].positions
    facecolor = 'green'
    if also_run_backwards and ix_frame > num_forward_frames:
        facecolor = 'purple'
    for i,circle in zip(xrange(init_container.num_particles), circles):
        circle.center = (posns[i][0], posns[i][1])  # x and y
        if i in special_particles:
            circle.set_facecolor('blue')
        else:
            circle.set_facecolor(facecolor)
    return circles

num_total_frames = num_forward_frames
if also_run_backwards:
    num_total_frames += num_forward_frames
frames = int(num_total_frames / frame_show_modulus)
anim = FuncAnimation(fig, next_frame, frames=frames, interval=dt, blit=True, init_func = init)
if save_animation:
    anim.save('pat_mol_dyn_{}_animation.avi'.format(sim_name), fps=30)
try:
    if show_animation:
        pl.show()
    else:
        pl.clf()
except:  # in true python style, ignore weird Tk error when closing plot window
    pass

# Plot PE
times = pl.frange(dt, num_total_frames*dt, dt)  # skip zeroth time because we have no value for it
pes = [c.potential_energy for c in containers[1:]]  # skip first container because it has no PE
plotted_pes = np.array(pes[:len(times)])
plotted_pes /= init_container.num_particles  # PE per particle
pl.plot(times, plotted_pes, 'o-', color='black', markersize=1, linewidth=0.1)
pl.ylabel('Potential energy per particle')
pl.xlabel('Time')
pl.title('PE/particle for {}'.format(sim_name))
pl.savefig('pat_mol_dyn_{}_pe.png'.format(sim_name))
pl.show()

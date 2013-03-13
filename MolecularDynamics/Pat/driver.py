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

from patku_sim import Container, VerletIntegrator
from problems import *

also_run_backwards = True
save_animation = True
num_forward_frames = 10 * 50
frame_show_modulus = 10  # only show every nth frame
dt = 1e-2
sim_name = ['symmetry', 'problem_3', 'problem_1', 'line'][0]
extra_params = {'Lx': 10,
    'random_velocity': 0.1, 'random_particle_ix': 0,
    'symmetry_number': 'six'}


# Circle code courtesy of Kevin Joyce
def get_nice_circle(x, y, radius = 0.5*particle_radius, color="lightsteelblue", facecolor="green", alpha=.6, ax=None ):
    """ add a circle to ax or current axes
    """
    e = pl.Circle([x, y], radius)
    if ax is None:
        ax = pl.gca()
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_edgecolor( color )
    e.set_linewidth(3)
    e.set_facecolor( facecolor )  # "none" not None
    e.set_alpha( alpha )
    return e

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
    e = get_nice_circle(posn[0], posn[1])
    circles.append(ax.add_patch(e))

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

frames = num_forward_frames
if also_run_backwards:
    frames = num_forward_frames*2 + 1  # include initial one
frames = int(frames / frame_show_modulus)
anim = FuncAnimation(fig, next_frame, frames=frames, interval=dt, blit=True)
if save_animation:
    anim.save('pat_mol_dyn_{}_animation.avi'.format(sim_name), fps=30)
pl.show()

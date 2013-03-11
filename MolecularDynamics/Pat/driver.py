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

particle_radius = 2**(1.0/6)

# Animation and circle code courtesy of Kevin Joyce
def circle(x, y, radius = 0.5*particle_radius, color="lightsteelblue", facecolor="green", alpha=.6, ax=None ):
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

init_container = Container()

##One particle, number 5, has a tiny velocity in the upwards direction and goes a tiny bit slower in x.
sim_name = 'line'
ix_particle_to_highlight = None
if sim_name == 'line':
    ix_particle_to_highlight = 5
    gamma = 1e-6
    Lx = 10
    Ly = 10
    init_container.bounds = (Lx, Ly)
    for i in xrange(1, Ly+1):
        position = [Lx / 2.0, (i-0.5)]
        if i == ix_particle_to_highlight:
            init_container.add_particle(position, [1-gamma, gamma])
        else:
            init_container.add_particle(position, [1, 0])


num_frames = 1000
dt = 1e-2
containers = [init_container]
integrator = VerletIntegrator()
for i in range(num_frames - 1):
    next_container = integrator.step(containers[-1], dt)
    containers.append(next_container)
end_container = containers[-1]

# Now run... backwards!
for i in range(num_frames - 1):
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

##particle_plots = [ax.plot([], [], marker='o')[0] for i in range(init_container.num_particles)]
## Pre initializing is necessary I guess
posns = init_container.positions
circles = []
for i,posn in enumerate(posns):
    facecolor = 'green'
    if i == ix_particle_to_highlight:
        facecolor = 'blue'
    e = circle(posn[0], posn[1], facecolor = facecolor)
    circles.append(ax.add_patch(e))

def next_frame(i):
##    global particle_plots
##    for i,plot in zip(xrange(init_container.num_particles), particle_plots):
##        plot.set_data(posns[i][0], posns[i][1])  # x and y
    posns = containers[i].positions
    for i,circle in zip(xrange(init_container.num_particles), circles):
        circle.center = (posns[i][0], posns[i][1])  # x and y
    return circles

# Twice the frames because we run backwards
anim = FuncAnimation(fig, next_frame, frames=num_frames*2 - 1, interval=dt, blit=True)
##anim.save('pat_mol_dyn_{}_animation.avi'.format(sim_name), fps=30)
pl.show()

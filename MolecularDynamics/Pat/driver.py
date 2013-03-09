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

from patku_sim import Container

container = Container(bounds=(5,5))
c = container
##One particle, number 5, has a tiny velocity in the upwards direction and goes a tiny bit slower in x.
initialization = 'line'
if initialization == 'line':
    gamma = 1e-6
    Lx = 10.
    Ly = 10.
    c.bounds = (Lx, Ly)
    for i in range(1, 12):
        position = [Lx / 2.0, (i-0.5)/11.0 * Ly]
        if i ==5:
            c.add_particle(position, [1-gamma, gamma])
        else:
            c.add_particle(position, [1, 0])
##container.accelerations = np.zeros((4,3))
##container.velocities = np.zeros((4,3))
##container.positions = np.array([
##    [],
##])
##container.add(list of particles, 2D or 3D)

# Animation code courtesy Kevin Joyce
########### FORCED ANIMATION ##########
########## PARAMS ###############
num_frames = 50
delay = 1
save_animation = True
file_name = "line"
#################################

circles = []
fig = pl.figure()
ax = pl.gca()
ax.set_aspect('equal')
ax.set_xlim((0, c.bounds[0]))
ax.set_ylim((0, c.bounds[1]))

def prettify_circle(e):
  color="lightsteelblue"
  facecolor="green"
  alpha=.6
  e.set_clip_box(ax.bbox)
  e.set_edgecolor( color )
  e.set_linewidth(3)
  e.set_facecolor( facecolor )  # "none" not None
  e.set_alpha( alpha )
  return e

## Pre initializing is necessary I guess
posns = c.positions
radius = 2**.2
for x in posns:
  e = Circle( (x[1], x[0]), radius=radius )
  e = prettify_circle(e)
  circles.append(ax.add_patch(e))

def init():
  return circles

def next_frame(i):
#  print "Frame: {}".format(i)
  for i in range(len(circles)):
    x = posns[i][0]
    y = posns[i][1]
    e = Circle((x, y), radius=radius)
    e.update_from(circles[i])
    circles[i] = ax.add_patch(e)
  return circles

anim = FuncAnimation(fig, next_frame, init_func=init, frames=num_frames, interval=delay, blit=True)

if save_animation:
  anim.save(file_name+".mp4",fps=25)
else:
  pl.show()


##import timeit
##
##timer = timeit.Timer(run)
##time_taken = timer.timeit(1)
##print 'Time', time_taken, 'seconds'
##
### Animate orbit
### Code courtesy of George Lesica
##fig = pl.figure(figsize=(8, 8))
##ax = pl.axes(xlim=(-2, 2), ylim=(-2, 2))  # necessary because initial plot is too zoomed in
##pl.title('{} Orbit Simulation'.format(name), fontsize=16)
##pl.xlabel('X Position')
##pl.ylabel('Y Position')
##p1, = ax.plot([], [], marker='o')
##p2, = ax.plot([], [], marker='o')
##p3, = ax.plot([], [], marker='o')
##
##def animate(i):
##    b1, b2, b3 = iter(systems[i])
##    p1.set_data([b1.x], [b1.y])
##    p2.set_data([b2.x], [b2.y])
##    p3.set_data([b3.x], [b3.y])
##    return p1, p2, p3
##
##anim = FuncAnimation(fig, animate, frames=len(times), interval=1, blit=True)
##anim.save('pat_{}_orbit_animation.avi'.format(name), fps=30)
####    pl.show()

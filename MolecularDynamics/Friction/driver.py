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

also_run_backwards = False
show_animation = True
save_animation = True
particle_radius = 2**(1.0/6)  # a.k.a. 'a'; actually, the point at which the Lenn-Jones force is stable between particles
frame_show_modulus = 10  # only show every nth frame
dt = 1e-2
sim_name = 'friction'
xlim, ylim = (0, 50), (-5, 5*particle_radius)


def create(cnt_sled_particles, cnt_floor_particles=100):
    """
    """
    global particle_radius
    c = Container()
    a = particle_radius
    # Make floor with particles spaced 'a' apart
    for i in xrange(cnt_floor_particles):
        c.add_particle([a*i, 0])
    # Make sled with particles in equilateral triangle of side length '2a' with two particles placed 'a' above floor starting at a/2
    # Connect sled particles like a trianglular truss
    return c

init_container = create(13, 5)
containers = [init_container]
integrator = VerletIntegrator()

# Animate orbit
# Code courtesy of George Lesica
fig = pl.figure(figsize=(8, 4))
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
while num_forward_frames < 500:
    num_forward_frames += 1
    for _ in xrange(1 + frame_show_modulus):  # always run at least once
        next_container = integrator.step(containers[-1], dt)
        containers.append(next_container)

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

anim = FuncAnimation(fig, next_frame, interval=dt, blit=True)
try:
    if show_animation:
        pl.show()
except:  # in true python style, ignore weird Tk error when closing plot window
    pass
pl.clf()

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
    anim.save('mol_dyn_friction_{}_animation.avi'.format(num_forward_frames), fps=30)


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

#-------------------------------------------------------------------------------
# Name:        orbits driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
import pylab as pl
import math
##import unittest
import itertools
from pprint import pprint, pformat
from collections import defaultdict, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+

from patku_sim import *


def run(file_path):
    global times, systems, integr
    # Driver idea
    system_init = FileReader.read_system_from(file_path)
    integr = Integrator(system_init.masses)
    t_start = 0
    t_end = 30
    dt = 0.1
    times = [t_start]
    systems = [system_init]
    while times[-1] < t_end:
        prev_time = times[-1]
        next_time = prev_time + dt
        next_system = integr.step(systems[-1], next_time, prev_time)
        times.append(next_time)
        systems.append(next_system)

import timeit
names_and_files = [
    ('Euler', r'C:\Users\Pat\Documents\My Dropbox\Simulation Hecuba Group Share\ObjectOrientedOrbits\euler_problem.txt'),
    ('Lagrange', r'C:\Users\Pat\Documents\My Dropbox\Simulation Hecuba Group Share\ObjectOrientedOrbits\lagrange_problem.txt'),
    ('Montgomery', r'C:\Users\Pat\Documents\My Dropbox\Simulation Hecuba Group Share\ObjectOrientedOrbits\montgomery_problem.txt')
]

for name, file_path in names_and_files:
    print 'Running simulation on', name, 'data'
    timer = timeit.Timer(lambda: run(file_path))
    time_taken = timer.timeit(1)
    print 'Time', time_taken, 'seconds'

    # Animate orbit
    # Code courtesy of George Lesica
    fig = pl.figure(figsize=(8, 8))
    ax = pl.axes(xlim=(-2, 2), ylim=(-2, 2))  # necessary because initial plot is too zoomed in
    pl.title('{} Orbit Simulation'.format(name), fontsize=16)
    pl.xlabel('X Position')
    pl.ylabel('Y Position')
    p1, = ax.plot([], [], marker='o')
    p2, = ax.plot([], [], marker='o')
    p3, = ax.plot([], [], marker='o')

    def animate(i):
        b1, b2, b3 = iter(systems[i])
        p1.set_data([b1.x], [b1.y])
        p2.set_data([b2.x], [b2.y])
        p3.set_data([b3.x], [b3.y])
        return p1, p2, p3

    anim = FuncAnimation(fig, animate, frames=len(times), interval=1, blit=True)
    anim.save('pat_{}_orbit_animation.avi'.format(name), fps=30)
##    pl.show()

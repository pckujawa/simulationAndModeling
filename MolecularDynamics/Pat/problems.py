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
##from matplotlib.animation import FuncAnimation  # v1.1+

from patku_sim import Container

# _extra_params = None
particle_radius = 2**(1.0/6)


def line_sim(container, special_particles, **extra_params):
    ##One particle, number 5, has a tiny velocity in the upwards direction and goes a tiny bit slower in x.
    ix_particle_to_highlight = 5
    special_particles.append(ix_particle_to_highlight)
    gamma = 1e-6
    Lx = 10
    Ly = 10
    container.bounds = (Lx, Ly)
    for i in xrange(1, Ly+1):
        position = [Lx / 2.0, (i-0.5)]
        if i == ix_particle_to_highlight:
            container.add_particle(position, [1-gamma, gamma])
        else:
            container.add_particle(position, [1, 0])
    # Maybe add assertions for distances between particles
    return container, special_particles

# Symmetry
# Try the following highly symmetric initial conditions to test your simulation. Two and four are probably the best tests, because it is easy to anticipate their behavior, but the others are interesting.
def symmetry(container, special_particles, **extra_params):
    p = container
    # distance and velocity here work well for Lx = 10.
##    Lx = extra_params.get('Lx', None)
##    Lx = Lx or container.bounds[0]  # error here would indicate bounds haven't been set
    Lx = 10
    container.bounds = (Lx, Lx)  # ?? What for Ly?
    dist = Lx / 5.
    vel = dist / 5.
    x,y = 5,5  # Jesse's stuff is centered at 5,5, whereas mine is at 0,0
    initialization = extra_params['symmetry_number']
    if initialization == 'one':
        p.add_particle([x, y+dist], [0, 0])
    elif initialization == 'two':
        p.add_particle([x - dist, y], [vel, 0.])
        p.add_particle([x + dist, y], [-vel, 0.])
    elif initialization == 'three':
        p.add_particle([x, y + dist*math.sqrt(3)/2], [0., -vel])
        p.add_particle([x - dist, y], [vel, 0.])
        p.add_particle([x + dist, y], [-vel, 0.])
    elif initialization == 'four':
        p.add_particle([x - dist, y], [vel, 0.])
        p.add_particle([x + dist, y], [-vel, 0.])
        p.add_particle([x, y + dist], [0., -vel])
        p.add_particle([x, y - dist], [0., vel])
    elif initialization == 'six':
        p.add_particle([x, y + dist], [0., -vel])
        p.add_particle([x, y - dist], [0., vel])
        other_x = x - dist/math.sqrt(2)
        x += dist/math.sqrt(2)
        p.add_particle([x, x], [-vel/math.sqrt(2), -vel/math.sqrt(2)])
        p.add_particle([other_x, x], [vel/math.sqrt(2), -vel/math.sqrt(2)])
        p.add_particle([other_x, other_x], [vel/math.sqrt(2), vel/math.sqrt(2)])
        p.add_particle([x, other_x], [-vel/math.sqrt(2), vel/math.sqrt(2)])
    elif initialization == 'eight':
        p.add_particle([x - dist, y], [vel, 0.])
        p.add_particle([x + dist, y], [-vel, 0.])
        p.add_particle([x, y + dist], [0., -vel])
        p.add_particle([x, x - dist], [0., vel])
        other_x = x - dist/math.sqrt(2)
        x += dist/math.sqrt(2)
        p.add_particle([x, x], [-vel/math.sqrt(2), -vel/math.sqrt(2)])
        p.add_particle([other_x, x], [vel/math.sqrt(2), -vel/math.sqrt(2)])
        p.add_particle([other_x, other_x], [vel/math.sqrt(2), vel/math.sqrt(2)])
        p.add_particle([x, other_x], [-vel/math.sqrt(2), vel/math.sqrt(2)])
    return container, special_particles


def problem_1(container, special_particles, **extra_params):
    """Create a triangular (hex) lattice with N = 64 particles. The dimensions of the lattice should be Lx and Ly = sqrt(3)/2*Lx
    """
    N = 8  # this is 8x8 = 64 particles
    Lx = extra_params.get('Lx', None)
    Lx = Lx or container.bounds[0]  # error here would indicate bounds haven't been set
    Ly = math.sqrt(3)/2 * Lx
    container.bounds = (Lx, Ly)
    dx = Lx / float(N)
    dy = Ly / float(N)
    add_triangle_lattice(container, N, dx, dy)
    return container, special_particles


def add_triangle_lattice(container, N, dx, dy):
    for i in range(N):
        for j in range(N):
            y = dy * (j + 0.5)
            if j % 2 == 0:
                x = dx * (i + 0.25)
            else:
                x = dx * (i + 0.75)
            container.add_particle([x, y], [0, 0])  # , mass=1


def problem_3(container, special_particles, **extra_params):
    """Initial positions of the particles are on the nodes of a square lattice. Choose N = 64 atoms and Lx = Ly = 9. All velocities are initially zero.
    """

    random_velocity = extra_params.get('random_velocity', None)
    random_particle_ix = extra_params.get('random_particle_ix', None)  # TODO `or random_ix`
    if random_particle_ix is not None:
        special_particles.append(random_particle_ix)
    N = 8  # this is 8x8 = 64 particles
    Lx = math.sqrt(math.sqrt(3)/2) * 9
    Ly = Lx
    container.bounds = (Lx, Ly)
    dx = Lx / float(N)
    dy = Ly / float(N)
    add_square_lattice(container, N, dx, dy, random_particle_ix, random_velocity)
    return container, special_particles

def add_square_lattice(container, N, dx, dy, random_particle_ix=None, random_velocity=None):
    for i in xrange(N):
        for j in xrange(N):
            vx, vy = 0, 0
            if i*8 + j == random_particle_ix:
                vy = random_velocity or 1  # TODO make random value
            x = dx * (i + 0.5)
            y = dy * (j + 0.5)
            container.add_particle([x, y], [vx, vy])  # , mass)


name_to_sim = {
    'line': line_sim,
    'problem_1': problem_1,
    'problem_3': problem_3,
    'symmetry': symmetry}


def get_container_for(s_sim_name, **extra_params):
    """Given the name of a simulation, return a container initialized to that sim's specs.
    :param extra_params: optional, named parameters to manipulate sims that allow it
    :returns: initialized container, list of special particles (e.g. to show in a different color)
    """
    global _extra_params
    # _extra_params = extra_params
    container = Container()
    special_particles = []
    return name_to_sim[s_sim_name](container, special_particles, **extra_params)

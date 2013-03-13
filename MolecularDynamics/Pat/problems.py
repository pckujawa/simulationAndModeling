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

_extra_params = None

def line_sim(container, special_particles):
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
    return container, special_particles

name_to_sim = {
    'line': line_sim }

def get_container_for(s_sim_name, **extra_params):
    """Given the name of a simulation, return a container initialized to that sim's specs.
    :param extra_params: optional, named parameters to manipulate sims that allow it
    :returns: initialized container, list of special particles (e.g. to show in a different color)
    """
    global _extra_params
    _extra_params = extra_params
    container = Container()
    special_particles = []
    return name_to_sim[s_sim_name](container, special_particles)

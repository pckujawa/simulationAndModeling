#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction driver
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
# import matplotlib as pl
import math
# from scipy.spatial import KDTree
##import unittest
##import itertools
##from pprint import pprint, pformat
# from collections import defaultdict#, namedtuple
#from matplotlib.animation import FuncAnimation  # v1.1+
#import timeit
#import cPickle as pickle

from patku_sim import Container, moldyn#, graphical, VerletIntegrator,
# from patku_sim.libs import Struct


def hourglass(funnel_width, hole_width, d_angle_from_horizon_to_wall, dist_between_anchors, **kwargs_ignore):
    """ Make a container with "anchor" particles set up as walls of an hourglass
    :param: kwargs_ignore ignored extra args so you can use an exploded dictionary (**myargdict) for arg input
    :returns: container, y_funnel_bottom
    """
    tan_over_2 = math.tan(math.radians(d_angle_from_horizon_to_wall))/2.0
    height_funnel = funnel_width * tan_over_2
    height_hole = hole_width * tan_over_2
    # Centers of wall particles
    def get_anchor_centers(xstart, xend):
        # np.arange doesn't let you go from a greater to a lesser value, so invert twice to get that behavior
        if xstart > xend:
            xs = -np.arange(-xstart, -xend, dist_between_anchors)
        else:
            xs = np.arange(xstart, xend, dist_between_anchors)
        return xs
    # Right wall equation: goes from ( (wf + wh)/2, hh - hf ) to (wf, 0)
    #m = (height_funnel - height_hole) / (funnel_width - (funnel_width + hole_width)/2.0)
    y_funnel_bottom = height_hole - height_funnel  # y is negative
    xrs = get_anchor_centers((funnel_width + hole_width)/2.0, funnel_width)
    yrs = get_anchor_centers(y_funnel_bottom, 0)
    # Left wall goes from (0, 0) to ( (wf - wh)/2, hh - hf)
    xls = get_anchor_centers(0, (funnel_width - hole_width)/2.0)
    yls = get_anchor_centers(0, y_funnel_bottom)

    # Place anchors
    c = Container()
    def add_anchors(xs, ys):
        for x, y in zip(xs, ys):
            c.add_particle([x, y])
    add_anchors(xls, yls)
    add_anchors(xrs, yrs)

    # Boundary so particles keep recirculating
    #  x is just a fudge value > funnel width so it doesn't affect anything
    #  y is some arbitrary line below the bottom of the funnel
##    c.bounds = (funnel_width+1, y_funnel_bottom - 5*hole_width)
    #TODO check that bounds code will allow negative y
    return c, y_funnel_bottom

def add_sand(container, funnel_width, dx, dy):
    moldyn.add_triangle_lattice(container, funnel_width, dx, dy)
    return container


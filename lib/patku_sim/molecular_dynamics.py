#-------------------------------------------------------------------------------
# Name:        molecular dynamics library
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
##import pylab as pl
import math
import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict, namedtuple

##from libs import Struct

##Point2D = namedtuple('Point2D', 'x y')
##Point3D = namedtuple('Point3D', 'x y z')

class Container(object):
    def __init__(self, bounds=None):
        self.bounds = bounds

    def bound(self, points):
        """Wrap points in space within this torus, ensuring that no dimension is out of bounds.
        :type points: numpy array
        :returns: points within bounds
        :rtype: numpy array
        """
        if not self.bounds:
            return points

        bounded = np.copy(points)
        # Because numpy doesn't handle multi-dimensional arrays the same as 1-dimensional ones, it's easiest to just make it always look like a multi-dim array
        cPoints, = points.shape  # can't unpack second item in 1-dim array
        if cPoints == 1:
            bounded = np.array([bounded, np.zeros_like(bounded)])
        ignore, cDims = bounded.shape
        for i in xrange(cDims):
            xs = bounded[:,i]
            boundary = self.bounds[i]
            xs[xs > boundary] -= boundary
            xs[xs < 0] += boundary
        if cPoints == 1:
            bounded = bounded[0]  # pull back out the 1-dim array
        return bounded


    def get_distance_matrices(self, *points):
        """Calculates the distances between all points for each dimension and radially.

        :type points: lists or arrays
        :returns: each dimension's distance matrix, in same order as passed in, followed by radial distance matrix
        """
        cPoints = len(points)
        if cPoints < 2:
            raise ValueError("Distance mtx for one point is the point's dimensions. Perhaps you meant to provide more than one point to this function. Maybe you need to unpack your list/tuple.")
        # Ensure each point has the same dimension
        cDim = len(points[0])  # count of dimensions
        for p in points:
            assert len(p) == cDim
        aPoints = np.array(points)
        # Use an inner iteration function because it's more versatile and easier to code than appending to a list
        def _iter():
            for i in xrange(cDim):
                xs = aPoints[:,i]
                xdist = np.tile(xs, (cPoints, 1))
                if self.bounds:
                    xbound = self.bounds[i]
                    xdist[xdist > xbound / 2.0] -= xbound
                    xdist[xdist < -xbound / 2.0] += xbound
                xdist = xdist - xdist.T
                yield xdist
        linear_distances = list(_iter())
        radial_distance = np.zeros_like(linear_distances[0])
        for x in linear_distances:
            radial_distance += x**2
        radial_distance = np.sqrt(radial_distance)
        linear_distances.append(radial_distance)  # too lazy to name a temp variable
        return linear_distances


class ContainerTests(unittest.TestCase):
    def test_given_body_at_origin_and_body_at_111_in_3D_distance_is_euclidean_for_each_dimension(self):
        b_origin = [0,0,0]
        b_111 = [1,1,1]
        # Expect mtxs for x, y, and z dimensions
        e_x = np.array([
            [0, 1],
            [-1, 0]])
        e_y = np.copy(e_x)
        e_z = np.copy(e_x)
        e_r = np.copy(e_x)*math.sqrt(3)  # sqrt of squares for radial dist
        f = Container(bounds=None)
        a_x, a_y, a_z, a_r = f.get_distance_matrices(b_origin, b_111)
        for e,a, in ((e_x, a_x), (e_y, a_y), (e_z, a_z)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_given_three_2D_points_all_distances_are_euclidean(self):
        p1 = [0,0]
        p2 = [1,1]
        p3 = [1,2]
        ex = np.array([
            [0, 1, 1],
            [-1, 0, 0],
            [-1, 0, 0]])
        ey = np.array([
            [0, 1, 2],
            [-1, 0, 1],
            [-2, -1, 0]])
        f = Container(bounds=None)
        ax, ay, ar = f.get_distance_matrices(p1, p2, p3)
        for e,a, in ((ex, ax), (ey, ay)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_given_2D_points_at_00_and_31_and_bounds_at_22_distances_are_11(self):
        ex = np.array([
            [0, 1],
            [-1, 0]])
        ey = np.copy(ex)
        f = Container(bounds=(2, 2))
        ax, ay, a_r = f.get_distance_matrices([0,0], [3,1])
        for e,a, in ((ex, ax), (ey, ay)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_two_points_near_boundaries_so_distance_wraps_around(self):
        f = Container(bounds=(1,))
        expected = np.array([
            [0, -0.2],
            [0.2, 0]])
        actual, a_r = f.get_distance_matrices([0.1], [0.9])
        for e,a in zip(expected.flat, actual.flat):
            self.assertAlmostEqual(e, a)

    #NOTE: A 'p' preceding a number in a test method name means ., e.g. 0p1 = 0.1
    def test_point_at_1p1_and_bounds_at_1_wraps_to_point_at_0p1(self):
        c = Container(bounds=(1,))
        expected = np.array([0.1])
        actual = c.bound(np.array([1.1]))
        self.assertAlmostEqual(expected[0], actual[0])



if __name__ == '__main__':
    unittest.main(verbosity=2)

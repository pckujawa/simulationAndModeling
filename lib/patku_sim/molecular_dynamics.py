#-------------------------------------------------------------------------------
# Name:        molecular dynamics library
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
import pylab as pl
import math
import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict, namedtuple

##from libs import Struct


# Circle code courtesy of Kevin Joyce
def get_nice_circle(x, y, radius, color="lightsteelblue", facecolor="green", alpha=.6, ax=None ):
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


def add_square_lattice(container, dimension, dx, dy, random_particle_ix=None, random_velocity=None):
    """
    :param dimension: length of one of two identical dimensions
    """
    for i in xrange(dimension):
        for j in xrange(dimension):
            vx, vy = 0, 0
            if i*8 + j == random_particle_ix:
                if random_velocity is not None:
                    vy = random_velocity
            x = dx * (i + 0.5)
            y = dy * (j + 0.5)
            container.add_particle([x, y], [vx, vy])  # , mass)


def add_triangle_lattice(container, dimension, dx, dy):
    """
    :param dimension: length of one of two identical dimensions
    """
    for i in xrange(dimension):
        for j in xrange(dimension):
            y = dy * (j + 0.5)
            if j % 2 == 0:
                x = dx * (i + 0.25)
            else:
                x = dx * (i + 0.75)
            container.add_particle([x, y], [0, 0])  # , mass=1


def lennardJonesForce(distance_matrices, radius=1.0, epsilon=1.0):
    """
    :param radius: same as $\sigma$
    :param epsilon: Constant of proportionality
    :returns: accelerations and potential energy
    :rtype: tuple(ndarray, float)
    """
    # Find position differences
    dr = distance_matrices[-1]

    # Eliminate zeros on diagonal
    dr[np.diag_indices_from(dr)] = 1

    # Build term to make it easy to get potential energy AND force:
    if np.any(dr == 0):
        raise RuntimeError("The distance between some particles is zero. That's bad.")
    over6 = (radius / dr)**6  # WARNING div by zero on diagonals yields inf
    over12 = over6**2

    # Force magnitude, from equation 8.3 of Gould
    magnitude = 24.0 * epsilon/dr * (2.0*over12 - over6)

    # Project onto components, sum all forces on each particle The dx,dy,dz
    # are zero on diagonal and save from dr = -1 on diagonal contributing to
    # force
    accelerations = []
    for dx in distance_matrices[:-1]:
        ax = np.sum(-magnitude * dx/dr, axis=1)
        accelerations.append(ax)
    accelerations = np.array(accelerations)

    # Easiest to compute energy at this point from Equation 8.2 of Gould
    # Note the triu to avoid double counting entries and the k=1 to avoid
    # self-potential due to dr being -1 to avoid divergence
    pe = 4.0 * epsilon * np.sum(np.sum(np.triu(over12 - over6, k=1)))
    return accelerations.transpose(), pe


class SledForcer(object):
    def __init__(self, equilibrium_distance, u, k = 10):
        """
        :param u: initial position of the bottom right sled particle
        """
        self.spring_eq_dist = equilibrium_distance
        self.k = k
        self.u = np.array(u)

    def apply_force(self, positions, time, damping_vs):
        """
        :param positions: list of lists/arrays, sled particle positions matching the convention that 0 is the bottom left and last is the bottom right
        :returns: accelerations for each particle in the sled due to springs, pulling, and damping
        :rtype: tuple(spring, pull, damp)
        """
        vx, vy = damping_vs
        size = len(positions)
        dms = get_distance_matrices(positions, one_point_ok=True)
        dr = dms[-1]
        force_spring = self.k * (dr - self.spring_eq_dist)  # NOTE: *Not* negative, as in Jesse's code!
        # Each sled particle is only connected to the particles within two indices of itself, so eliminate others (and self-effects)
        #  There's probably a slicker way to create a matrix for the connections...
        if size > 1:
            diag1 = np.ones((size-1,))
            diag2 = np.ones((size-2,))
            spring_connections = \
                np.diagflat(diag1, 1) + \
                np.diagflat(diag1, -1) + \
                np.diagflat(diag2, 2) + \
                np.diagflat(diag2, -2)
            no_spring_connections = spring_connections != 1
            force_spring[no_spring_connections] = 0
        # Finally, find the directional accelerations
        spring_accels = np.zeros_like(positions)
        damping_as = np.zeros_like(positions)  # only use the first slot
        pulling_accels = np.zeros_like(positions)  # only use the last slot
        dr += np.eye(*dr.shape)  # 'hide' zeros on diagonal from denominator
        displacement = positions[-1] - self.u
        for ix_dim, dx in enumerate(dms[:-1]):  # skip dr
            x_hat = dx/dr  # not exactly x^, but similar
            directional_accels = x_hat * force_spring
            spring_accels[:, ix_dim] = np.sum(directional_accels, axis=1)  # one accel per particle per dimension
            damping_as[0, ix_dim] = damping_vs[ix_dim]
            if ix_dim == 0:  # only horizontal pulling allowed
                pulling_accels[-1, ix_dim] = (self.pulling_force_v*time - displacement[ix_dim])
        damping_as *= -self.damp_force_multiplier
        if not self.allow_negative_pull_force:
            pulling_accels[pulling_accels < 0] = 0
        pulling_accels *= self.pulling_force_k
        # Normal force
        num_particles = len(positions)
        normal_as = np.zeros_like(positions)
        normal_as[:, 1] = -self.normal_force/num_particles  # only y component
        return spring_accels, pulling_accels, damping_as, normal_as


class Container(object):
    def __init__(self, bounds=None):
        self.bounds = bounds
        self._positions = []
        self._velocities = []
        self._accelerations = []
        self.kinetic_energies = []
        self.potential_energy = None
        self.num_particles = 0
        self.time = 0  # point in time

    def apply_force(self, sled_forcer=None):
        distance_matrices = self.get_distance_matrices()
        accelerations, pe = lennardJonesForce(distance_matrices)
        self.potential_energy = pe
        try:
            # If we have a sled, apply spring forces
            ix_damper = self.sled_particle_ixs[0]
            ix_puller = self.sled_particle_ixs[-1]
            sled_posns = self._positions[ix_damper:]
            assert ix_puller == self.num_particles-1
            damping_vs = self._velocities[ix_damper]  # first sled particle
            spring, pull, damp, normal = sled_forcer.apply_force(sled_posns, self.time, damping_vs)
            self.spring_accelerations = spring
            self.pull_accelerations = pull[-1]  # one particle, but all dims
            self.damp_accelerations = damp[0]  # same
            sled_accelerations = spring + pull + damp + normal
            accelerations[self.sled_particle_ixs] += sled_accelerations
            # If we have a floor, keep it still
            accelerations[self.floor_particle_ixs] = 0  # broadcasts across all dimensions
        except AttributeError:
            pass
        self._accelerations = accelerations.tolist()

    def add_particle(self, position, velocity=None, acceleration=None):
        """Add one particle of n dimensions, *unbounded*.
        """
        self._positions.append(np.array(position))
        v = velocity or np.zeros_like(position)
        a = acceleration or np.zeros_like(position)
        self._velocities.append(np.array(v))
        self._accelerations.append(np.array(a))
        self.num_particles += 1

    @property
    def positions(self):
        return np.array(self._positions)

    @property
    def velocities(self):
        return np.array(self._velocities)

    @property
    def accelerations(self):
        return np.array(self._accelerations)

    def set_velocities(self, velocities):
        """
        :type velocities: ndarray or list
        """
        if type(velocities) is list:
            self._velocities = velocities
        else:
            self._velocities = velocities.tolist()

    def set_positions_velocities_accelerations(self, positions, velocities=None, accelerations=None):
        """Bound points
        :type positions: ndarray
        :type velocities: optional ndarray
        :type accelerations: optional ndarray
        """
        shape = positions.shape  # stores num dimensions and num particles
        self.num_particles = shape[0]
        assert self.num_particles == len(positions)
        self._positions = self.bound(positions).tolist()
        if velocities is None:
            velocities = np.zeros_like(positions)
        self.set_velocities(velocities)
        if accelerations is None:
            accelerations = np.zeros_like(positions)
        self._accelerations = accelerations.tolist()

    def bound(self, points):
        """Wrap points in space within this torus, ensuring that no dimension is out of bounds.
        :type points: numpy array
        :returns: points within bounds
        :rtype: numpy array
        """
        if self.bounds is None:
            return points

        bounded = np.copy(points)
        # Because numpy doesn't handle multi-dimensional arrays the same as 1-dimensional ones, it's easiest to just make it always look like a multi-dim array
        points_shape = points.shape
        cPoints = points_shape[0]
        if cPoints == 1:
            bounded = np.array([bounded, np.zeros_like(bounded)])
        _ignore, cDims = bounded.shape
        for i in xrange(cDims):
            xs = bounded[:,i]
            boundary = self.bounds[i]
            # TODO place in while loop to ensure positions haven't gone multiple times past the boundary
            xs = xs % boundary
##            xs[xs > boundary] -= boundary
##            xs[xs < 0] += boundary
            bounded[:,i] = xs  # is this necessary? seems so
        if cPoints == 1:
            bounded = bounded[0]  # pull back out the 1-dim array
        return bounded

    def get_distance_matrices(self):
        return get_distance_matrices(self._positions, self.bounds)


def get_distance_matrices(points, bounds=None, one_point_ok=False):
    """Calculates the distances between all points for each dimension and radially.
    :param bounds: tuple, linear boundaries for each dimension, assuming origin at 0
    :type points: list or ndarray
    :returns: each dimension's distance matrix, in same order as passed in, followed by radial distance matrix
    """
    cPoints = len(points)
    if cPoints < 2 and not one_point_ok:
        raise ValueError("Distance mtx for one point is the point's dimensions. Perhaps you meant to provide more than one point to this function. Maybe you need to unpack your list/tuple.")
    # Ensure each point has the same dimension
    cDim = len(points[0])  # count of dimensions
    for p in points:
        assert len(p) == cDim
    aPoints = np.array(points)
    # Use an inner iteration function because it's more versatile and easier to code than appending to a list
    def _iter():
        for i in xrange(cDim):
            xs = aPoints[:, i]
            xdist = np.tile(xs, (cPoints, 1))
            xdist = xdist - xdist.T
            if bounds is not None:
                try:
                    xbound = bounds[i]
                except IndexError:
                    raise Exception("There aren't enough boundaries for the number of dimensions in the points. Ensure that your bounds are of the same dimension as your points.")
                xdist[xdist > xbound / 2.0] -= xbound
                assert not np.any(xdist > xbound/2.0)
                xdist[xdist < -xbound / 2.0] += xbound
                assert not np.any(xdist < -xbound/2.0)
            yield xdist
    linear_distances = list(_iter())
    radial_distance = np.zeros_like(linear_distances[0])
    for x in linear_distances:
        radial_distance += x**2
    radial_distance = np.sqrt(radial_distance)
    linear_distances.append(radial_distance)  # too lazy to name a temp variable
    return linear_distances


class VerletIntegrator(object):
    def step(self, container, dt):
        """Move the container's particles forward in time and encapsulate in a new container.
        """
        c = container
        c_velocities = c.velocities
        c_accelerations = c.accelerations
        c_next = Container(bounds=c.bounds)
        # KLUDGE copying attributes should be handled in Container to prevent inconsistencies
        c_next.time = c.time + dt
        c_next.num_particles = c.num_particles
        try:
            c_next.floor_particle_ixs = c.floor_particle_ixs
            c_next.sled_particle_ixs = c.sled_particle_ixs
        except AttributeError:
            pass
        def _posn_iter():
            for xs, vxs, axs in zip(c.positions, c.velocities, c_accelerations):
                yield xs + vxs*dt + 0.5*axs*dt**2
        new_posns = list(_posn_iter())

        # Need to hold on to previous velocities
        c_next.set_positions_velocities_accelerations(
            np.array(new_posns), c_velocities)
        try:
            c_next.apply_force(self.sled_forcer)
        except AttributeError:
            c_next.apply_force()
        kes = None  # Kinetic Energy
        lVxs = []  # l=list; store vxs values
        for xs, vxs, axs, axs_prev in zip(
                c_next.positions, c_next.velocities, c_next.accelerations, \
                c_accelerations):
            new_vxs = vxs + 0.5*(axs + axs_prev)*dt
            lVxs.append(new_vxs)
            if kes is None:
                kes = np.zeros_like(new_vxs)
            kes += 0.5 * new_vxs**2  # mass is 1
        c_next.set_velocities(np.array(lVxs))
        c_next.kinetic_energies = kes
        return c_next


class VerletIntegratorTests(unittest.TestCase):
##positions should change due to both velocity and acceleration
##container result should have same bounds as initial container
##container result should have positions bound
##ideally, sum of PE and KE is constant throughout
    pass



class DistanceMatrixTests(unittest.TestCase):
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
        a_x, a_y, a_z, a_r = get_distance_matrices([b_origin, b_111])
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
        ax, ay, ar = get_distance_matrices([p1, p2, p3])
        for e,a, in ((ex, ax), (ey, ay)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_given_2D_points_at_00_and_31_and_bounds_at_22_distances_are_11(self):
        ex = np.array([
            [0, 1],
            [-1, 0]])
        ey = np.copy(ex)
        ax, ay, a_r = get_distance_matrices([[0,0], [3,1]], bounds=(2, 2))
        for e,a, in ((ex, ax), (ey, ay)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_two_points_near_boundaries_so_distance_wraps_around(self):
        expected = np.array([
            [0, -0.2],
            [0.2, 0]])
        actual, a_r = get_distance_matrices([[0.1], [0.9]], bounds=(1,))
        for e,a in zip(expected.flat, actual.flat):
            self.assertAlmostEqual(e, a)

class ContainerTests(unittest.TestCase):
    #NOTE: A 'p' preceding a number in a test method name means ., e.g. 0p1 = 0.1
    def test_point_at_1p1_and_bounds_at_1_wraps_to_point_at_0p1(self):
        c = Container(bounds=(1,))
        expected = np.array([0.1])
        actual = c.bound(np.array([1.1]))
        self.assertAlmostEqual(expected[0], actual[0])



if __name__ == '__main__':
    unittest.main(verbosity=2)

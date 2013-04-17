#-------------------------------------------------------------------------------
# Name:        molecular dynamics library
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
import pylab as pl
from scipy.spatial import KDTree
import scipy.sparse
import math
import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict, namedtuple

from libs import Struct


class Force(Struct):
    name = 'Unimplemented'
    attrs_to_copy = []
    def __init__(self, **other_attributes):
##        self.attrs_to_copy = attrs_to_copy
        Struct.__init__(self, **other_attributes)

    def __call__(self, container, sim_wide_params):
        raise NotImplementedError('Subclasses of Force must implement __call__(self, container)')



def get_pairs_within_distance(container, dist):
    posns = container.positions
    kd = scipy.spatial.KDTree(posns)
    return kd.query_pairs(dist)


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


def add_triangle_lattice(container, xlim, ylim, dx, dy):
    """Add particles to a container at x,y > 0 in a triangular configuration.
    """
    xs_base = np.arange(xlim[0], xlim[1], dx) + 0.25
    xs_odd = xs_base + 0.5
    ys_base = np.zeros_like(xs_base)
    ysteps = np.arange(ylim[0], ylim[1], dy) + 0.5
    for ix, y in enumerate(ysteps):
        ys = ys_base + y
        xs = xs_base
        if ix % 2 == 0:
            xs = xs_odd  # use cached values
        for x, y in zip(xs, ys):
            container.add_particle([x, y], [0, 0])

def lj_given_pair(pair_posns, bounds=None):
    """
    :returns: acceleration (for each dim) felt by particles, in order
    """
##    f, s = pair_posns  # first, second (n-dim lists/array of positions)
    # Find distance between the two
    dms = get_distance_matrices(pair_posns, bounds)
    return lennardJonesForce(dms)

# From KevinJ
def sparse_tile(x_positions, sp_dx):
    """
    Make a sparse matrix like sp_dx but with data from x_positions.
    :type x_positions: CSC sparse matrix
    """
    assert scipy.sparse.isspmatrix_csc(sp_dx)
    sp_x = sp_dx.copy()
    for i,j in zip(sp_dx.indptr[:-1], sp_dx.indptr[1:]):
        sp_x.data[i:j] = x_positions[sp_dx.indices[i:j]]
    return sp_x

def lj_given_kd_tree(kd, positions, distance=2.5):
    raise NotImplementedError()
    # TODO CSC sparse matrices can't do most of the things we need to do in LJ (fancy indexing, division). Maybe another sparse matrix type can?
    c = containers[-1]
    positions = c.positions
    kd = KDTree(positions)
    sm_dr = kd.sparse_distance_matrix(kd, distance)
##    sm_dr.getrow(0)
##    kd.query(positions[0], k=10, distance_upper_bound=2.5)
    sm_dr_csc = sm_dr.tocsc()  # F's up fancy ixing and multiplying
    sm_dxs = []
    for dim_ix in xrange(positions.shape[1]):  # KLUDGE, fix with other code that iterates over dims
        smx = sparse_tile(positions[:,dim_ix], sm_dr_csc)
        smdx = smx - smx.transpose()
        sm_dxs.append(smdx.todok())
    distance_matrices = sm_dxs + [sm_dr]
    radius=1.0; epsilon=1.0
    return lennardJonesForce(distance_matrices)


def lennardJonesForce(distance_matrices,
        radius=1.0, epsilon=1.0, cutoff_dist=float('inf'),
        anchor_ixs = None):
    """
    :param radius: same as $\sigma$
    :param epsilon: Constant of proportionality
    :param cutoff_dist: distance beyond which to ignore this force
    :returns: accelerations and potential energy
    :rtype: tuple(ndarray, float)
    """
    # Find position differences
    dr = distance_matrices[-1]
    ixs_not_touching = dr > cutoff_dist

    if anchor_ixs is not None:
        # How can we make sure anchors' effects aren't felt on other anchors?
##        assert dr[:, anchor_ixs][anchor_ixs].shape == tuple([len(anchor_ixs)]*2)
        dr[:, anchor_ixs][anchor_ixs] = 1e10  # far enough away to have insignificant effect

    # Eliminate zeros on diagonal
    dr[np.diag_indices_from(dr)] = 1

    if np.any(dr == 0):
        raise RuntimeError("The distance between some particles is zero. That's bad.")

    # If we're using sparse matrices, the radius and epsilon can't be floats because division by spm isn't supported. So turn them into spms.
    if scipy.sparse.isspmatrix(dr):
##        from scikits.learn.preprocessing import Normalizer
        def sparse_like(spmatrix, value):
            x = spmatrix.copy()
            x[spmatrix.nonzero()] = value  # SLOW!!!
            return x
        radius = sparse_like(dr, radius)
        epsilon = sparse_like(dr, epsilon)

    # Build term to make it easy to get potential energy AND force:
    over6 = (radius / dr)**6
    over12 = over6**2

    # Force magnitude, from equation 8.3 of Gould
    if np.isinf(cutoff_dist):
        magnitude = 24.0 * epsilon/dr * (2.0*over12 - over6)
    else:
        magnitude = 24.0 * epsilon/dr * (2.0*over12)  # no attraction

    # Project onto components, sum all forces on each particle The dx,dy,dz
    # are zero on diagonal and save from dr = -1 on diagonal contributing to
    # force
    accelerations = []
    for dx in distance_matrices[:-1]:
        dx[ixs_not_touching] = 0
        ax = np.sum(-magnitude * dx/dr, axis=1)
        accelerations.append(ax)
    accelerations = np.array(accelerations)

    # Easiest to compute energy at this point from Equation 8.2 of Gould
    # Note the triu to avoid double counting entries and the k=1 to avoid
    # self-potential due to dr being -1 to avoid divergence
    diff_term = over12 - over6
    diff_term[ixs_not_touching] = 0
    pe = 4.0 * epsilon * np.sum(np.sum(np.triu(diff_term, k=1)))
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


#TODO Make a container Builder for adding particles so that pulling out positions etc is a one-time cost
class Container(object):
    def __init__(self, bounds=None):
        self.bounds = bounds
        self._positions = []
        self._velocities = []
        self._accelerations = []
        self.kinetic_energies = []
        self.potential_energy = None
        self.num_particles = 0
        self.time = 0  # point in time that this snapshot belongs to

    def apply_force(self, forces=None, sim_wide_params=None, neighbor_facilitator=None):
        accel_map = {}
        if neighbor_facilitator is not None:
            nf = neighbor_facilitator
            posns = self.positions
##            accelerations, pe = lj_given_kd_tree(nf.generate_kdtree_if_needed(posns), posns, nf.d)
            pairs = nf.generate_pairs_if_needed(posns)
            acc_as = np.zeros_like(posns)  # accumulator for accels
            pe = 0
            for pair in pairs:
                pair = list(pair)  # needed for indexing
                pair_posns = posns[pair]
                accels, per_pe = lj_given_pair(pair_posns)
                acc_as[pair] += accels
                pe += per_pe
            accelerations = acc_as
##            # Compare
##            n_as, n_pe = accelerations, pe
##            distance_matrices = self.get_distance_matrices()
##            dm_as, dm_pe = lennardJonesForce(distance_matrices)
##            diff_pe = dm_pe - n_pe
##            assert abs(diff_pe) < 1
##            assert np.allclose(n_as, dm_as, atol=0.1)
##            diff_as = dm_as - n_as
##            print "Diff a,pe = {},{}".format(diff_as, diff_pe)
        elif forces is not None:
            for force in forces:
                accelerations = force(self, sim_wide_params)
                accel_map[force.name] = accelerations
                # If the force specifies attributes that should be copied to this container (e.g. LJ would have PE), do so
                try:
                    for attr in force.attrs_to_copy:
                        setattr(self, attr, getattr(force, attr))
                except AttributeError:
                    pass
        else:
            distance_matrices = self.get_distance_matrices()
            accelerations, pe = lennardJonesForce(distance_matrices)
            self.potential_energy = pe
##            accelerations[self.sled_particle_ixs] += sled_accelerations
##            # If we have a floor, keep it still
##            accelerations[self.floor_particle_ixs] = 0  # broadcasts across all dimensions
##        for a in accel_map.itervalues():
##            if np.any(a > 1):
##                pass
        accelerations = np.sum(accel_map.values(), axis=0)  # sum dim values

        try:
            anchor_ixs = sim_wide_params.anchor_ixs
            # http://stackoverflow.com/questions/7741878/how-to-apply-numpy-linalg-norm-to-each-row-of-a-matrix
            self.anchor_accels = np.sum(np.abs(accelerations[anchor_ixs])**2, axis=-1)**(1./2)
            accelerations[anchor_ixs] = 0  # all dims
        except AttributeError: pass
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

        bounded = np.copy(points)  # SLOW!!!
        # Because numpy doesn't handle multi-dimensional arrays the same as 1-dimensional ones, it's easiest to just make it always look like a multi-dim array
        points_shape = points.shape
        cPoints = points_shape[0]
        if cPoints == 1:
            bounded = np.array([bounded, np.zeros_like(bounded)])
        _ignore, cDims = bounded.shape
        for i in xrange(cDims):
            xs = bounded[:,i]
            min_b, max_b = self.bounds[i]
            assert min_b < max_b
            width = max_b - min_b
            # (EDIT: Wrong!) Because of the way that mod works (it wraps negative values around, rather than returning -(abs(x) % abs(y))), we can just use it straight
            # Need to treat neg and pos values different because of behavior of mod operator
            # On second thought, don't use mod, just assume small jumps
            too_far_neg_ixs = xs < min_b
            xs[too_far_neg_ixs] += width
            too_far_pos_ixs = xs > max_b
            xs[too_far_pos_ixs] -= width
            bounded[:,i] = xs  # is this necessary? seems so
        if cPoints == 1:
            bounded = bounded[0]  # pull back out the 1-dim array
        return bounded

    def dims_iter(self):
        """ Dimensions iterator (numerical ixs)
        """
        return xrange(len(self._positions[0]))

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
                    min_b, max_b = bounds[i]
                    width = max_b - min_b
                except IndexError:
                    raise Exception("There aren't enough boundaries for the number of dimensions in the points. Ensure that your bounds are of the same dimension as your points.")
                # Can't use mod because the lower triangle is negative and it wraps around weird
##                xdist = xdist % (width / 2.0)
                xdist[xdist > width / 2.0] -= width
                assert not np.any(xdist > width/2.0)
                xdist[xdist < -width / 2.0] += width
                assert not np.any(xdist < -width/2.0)
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
            c_next.apply_force(self.forces, sim_wide_params = self.sim_wide_params)
        except AttributeError:
            try:
                c_next.apply_force(self.sled_forcer, self.neighbor_facilitator)
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
        ax, ay, a_r = get_distance_matrices([[0,0], [3,1]], bounds=[(0, 2),(0, 2)])
        for e,a, in ((ex, ax), (ey, ay)):
            self.assertTrue(np.all(np.equal(e, a)))

    def test_two_points_near_boundaries_so_distance_wraps_around(self):
        expected = np.array([
            [0, -0.2],
            [0.2, 0]])
        actual, a_r = get_distance_matrices([[0.1], [0.9]], bounds=[(0, 1)])
        for e,a in zip(expected.flat, actual.flat):
            self.assertAlmostEqual(e, a)

class ContainerTests(unittest.TestCase):
    #NOTE: A 'p' preceding a number in a test method name means ., e.g. 0p1 = 0.1
    def get_bound_value(self, container, scalar):
        return container.bound(np.array([scalar]))[0]

    def test_point_at_1p1_and_bounds_at_1_wraps_to_point_at_0p1(self):
        c = Container(bounds=[(0, 1)])
        expected = np.array([0.1])
        actual = c.bound(np.array([1.1]))
        self.assertAlmostEqual(expected[0], actual[0])

    def test_min_and_max_bounds(self):
        l, r = (-1, 2)
        c = Container(bounds=[(l, r)])
        # out of bounds neg
        self.assertAlmostEqual(r - l - 1.9, self.get_bound_value(c, -1.9))
        # within bounds pos
        self.assertAlmostEqual(1.1, self.get_bound_value(c, 1.1))
        # within bounds neg
        self.assertAlmostEqual(-0.9, self.get_bound_value(c, -0.9))


if __name__ == '__main__':
    unittest.main(verbosity=2)

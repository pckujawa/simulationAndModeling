#-------------------------------------------------------------------------------
# Name:        patku_sim library
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
import pylab as pl
import math
import unittest
import itertools
from pprint import pprint, pformat
from collections import defaultdict
##from matplotlib.animation import FuncAnimation  # v1.1+

# http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
def chunk_iter(l, n, strict=False):
    """ Yield successive n-sized chunks from l.
        NOTE: This function return an iterator, not a list.
        You can either put it in a for loop or call
        `list(chunk_iter(...))`.
    Example:
states = np.array([1,1,1,1, 2,2,2,2, 3,3,3,3])
for x, y, vx, vy in chunk_iter(states, 4):
    print x  # could also access y, vx, vy
    """
    for i in xrange(0, len(l), n):
        if strict and i+n > len(l):
            raise ValueError("Not enough values to make a full chunk")
        yield l[i:i+n]


class Struct(object):
    '''An expando object.'''
    def __init__(self, **kwargs):
        '''Create a new object with attributes specified by the keyword args.
        Example:
            myobj = Struct(x=0, y='y', z=1.5)
            assert myobj.x == 0
            assert myobj.y == 'y'
            assert myobj.z == 1.5
        '''
        self.__dict__.update(kwargs)

    def __str__(self):
        '''Return a pretty string representation of this object.'''
        return str(self.__dict__)

    def __repr__(self):
        '''Return a pretty string representation of this object.'''
        return self.__str__()


class Body(Struct):
    '''stores, and retrieves the state data for one particle. In two dimensions, this consists of , postions, velocities, and mass respectively. Methods for getting and setting the state should be included. It will be helpful if the get_state function returns data in a way that lends itself to the conversion of state data to the sorts of vectors that are anticipated by the ode integration (see System, below).
    '''
    counter = itertools.count(0)  # for assigning unique ids

    def __init__(self, id_=None, **kwargs):
        super(type(self), self).__init__(**kwargs)
        if id_ is not None:
            self.id = id_
        else:
            self.id = Body.counter.next()


class System(Struct):
    '''This is a wrapper for a set of body objects. A list of bodies worked well for me. There is minimal other data, the total number of particles is needed. There should be a method for getting the state, and that method should return the data in a manner that is compatible with the vectors needed by the integration code. There should also be a way to set the state, given the vector output of the ode integrator. Finally, add functionality to add or remove a single particle.
    '''
    def __init__(self, bodies, **kwargs):
        '''Create a system with a copy of the arg bodies.'''
        self.bodies = bodies[:]  # copy
        self.bodies.sort(key = lambda b: b.id)  # maintain in sorted order
        self.num_bodies = len(self.bodies)
        super(type(self), self).__init__(**kwargs)

    def __iter__(self):
        return iter(self.bodies)

    @property
    def masses(self):
        return [b.mass for b in self.bodies]

    @masses.setter
    def masses(self, values):
        if len(values) != len(self.bodies):
            raise ValueError("Can't assign masses to the system's bodies because their lengths mismatch.")
        for b,v in zip(self.bodies, values):
            b.mass = v


##    def __setattr__(self, name, val):
##        # If input is a collection, try to broadcast to bodies
##        if type(val) in (list, np.array):
##            if len(val) != len(self.bodies):
##                raise ValueError("Can't broadcast the attribute '{}' to the system's bodies because their lengths mismatch.".format(name))
##            for body, item in zip(self.bodies, val):
##                setattr(body, name, item)
##        else:
##            return NotImplemented


class SystemStateAdapter(object):
    state_attributes = ['x', 'y', 'vx', 'vy']

    @classmethod
    def body_from_state(cls, state):
        b = Body()
        # Assign x,y, etc to the body, pulling from array
        for attr, val in zip(SystemStateAdapter.state_attributes, state):
            setattr(b, attr, val)
        return b

    @classmethod
    def state_from_body(cls, body):
        state = []
        for attr in SystemStateAdapter.state_attributes:
            val = getattr(body, attr)
            state.append(val)
        # Could be more concise but less debug-able:
##        return [getattr(body, attr) for attr in SystemStateAdapter.state_attributes]
        return state

    @classmethod
    def system_from_state(cls, state):
        """

        :param state: all bodies' positions and velocities as [x1,y1,vx1,vy1, x2,y2,vx2,vy2, ...]
        :type state: big 1D array or list
        """
        bodies = []
        for state_chunk in chunk_iter(state, len(SystemStateAdapter.state_attributes), True):
            b = SystemStateAdapter.body_from_state(state_chunk)
            bodies.append(b)
        return System(bodies)

    @classmethod
    def state_from_system(cls, system):
        cumulative_state = []
        for body in system:
            # NOTE: `extend`, not `append`
            cumulative_state.extend(SystemStateAdapter.state_from_body(body))
        return cumulative_state


class SystemStateAdapterTests(unittest.TestCase):
    def test_converts_state_to_individual_properties(self):
        e = Struct(x=1, y=2, vx=3, vy=4)  # expected
        state = np.array([e.x, e.y, e.vx, e.vy])
        target = SystemStateAdapter.body_from_state(state)
        for attr in SystemStateAdapter.state_attributes:
            self.assertEqual(getattr(e, attr),
                getattr(target, attr), attr + ' was different')



    # Initialize a holder for forces to later verify that, e.g. F[1][2] = -F[2][1]
_forces_table = {'x': defaultdict(dict), 'y': defaultdict(dict)}


class GravitationalForcer(object):
    """This class should provide right hand side functions  of the sort used by the solvers in scipy; With the normal cautions about f(t,x) vs. f(x,t) (for ode). This would be a good candidate for inheritance, implementing a basic interface to the , a conversion to go from  to , and a way of handling parameters. Derived classes would then implement the specific , and name the parameters.
    """
    def __init__(self, system, time_is_first_arg=True):
        self.system = system
        self.time_is_first_arg = time_is_first_arg

    def __call__(self, *args):
        '''RHS of the step function; provides the derivatives in order.
            :param args: either time followed by state (if self.time_is_first_arg is true) or state followed by time
            :type args: tuple or list
        '''
        if self.time_is_first_arg:
            return self._ode_fn(*args)  # unpack args in same order passed in
        return self._ode_fn(args[1], args[0])  # switch order of args

    def _ode_fn(self, t, state):
        # NOTE: Throughout this process, order must be maintained
        # Unpack state into bodies/system
        SystemBuilder().with_state(state).with_masses(masses).get()
        # Find acceleration of bodies
        # Return velocities and accels
        print 't:', t, 'state:', state


class GravitationalForcerTests(unittest.TestCase):
    def test_given_state_return_derivatives_in_order(self):
        target = GravitationalForcer(None)
        target(None, TODO)

def get_bodies_with_acceleration(bodies):
    G = 1  # not true, but fine for testing
    # NOTE: Instead of nested for loop, we can use itertools.permutation()
    for first in bodies:
        axs, ays = [], []  # all accelerations for first body
        for second in bodies:
            if first is second:
                continue  # don't look at force of body on itself
            m1, m2 = first.mass, second.mass
            x1, x2 = first.x, second.x
            y1, y2 = first.y, second.y
            r12x, r12y = x2 - x1, y2 - y1
            r12mag = np.sqrt(r12x**2 + r12y**2)
            # NOTE: accel *experienced* by first body
            get_accel = lambda r_component: G * m2 * r_component / r12mag**3
            a12x = get_accel(r12x)
            a12y = get_accel(r12y)
            axs.append(a12x)
            ays.append(a12y)
            _forces_table['x'][first.id][second.id] = a12x * m1
            _forces_table['y'][first.id][second.id] = a12y * m1
        first.ax, first.ay = sum(axs), sum(ays)
    return bodies  # mutated, so caller can instead just use her reference


class MiscTests(unittest.TestCase):
    def test_accelerations(self):
        b1 = Body(x=0, y=0, mass=1)
        b2 = Body(x=1, y=0, mass=2)
        target = System([b1, b2])
        get_bodies_with_acceleration(target.bodies)
        self.assertGreater(b1.ax, 0)  # pulled right
        self.assertAlmostEqual(b1.ax * b1.mass, -b2.ax * b2.mass)  # opposing forces

class Integrator(object):
    '''Work is done here. The integration machinery should be established, functions bound, and the forward integration carried out by methods in this class.
    '''
    pass

class IntegratorTests(unittest.TestCase):
    def test_next_returns_next_time_step(self):
        states = odeint(double_star, double_star_state_0, times)
        TODO


class FileReader(object):
    '''Provides factory methods for reading and converting files.'''
    @classmethod
    def read_system_from(cls, file_path):
        '''Read a file and convert contents into a System.
        File contents are of the form:
            N = Number of bodies
            x1 y1 vx1 vy1 m1
            x2 y2 vx2 vy2 m2
            ...

        Arguments:
        :param file_path: path to the file to read
        :type file_path: string
        :returns: System represented by file contents
        :rtype: System
        '''
        f = file(file_path)
        f.next()  # line with number of bodies - ignore
        lines = np.genfromtxt(f, delimiter=' ')
        # Pull out the state values (x,y,vx,vy) and combine them
        state = lines[:,0:4].flatten()
        masses = lines[:,4]
        s = SystemStateAdapter.system_from_state(state)
        s.masses = masses
        return s


def run():
    # Driver idea
    system_init = FileReader.read_system_from('path/to/file.txt')
    integr = Integrator()
    t_start = 0
    t_end = 10
    dt = 0.1
    times = [t_start]
    systems = [system_init]
    while times[-1] < t_end:
        next_time = times[-1] + dt
        force_fn = GravitationalForcer(system)
        next_system = integr.next(force_fn, systems[-1], next_time)
        times.append(next_time)
        systems.append(next_system)


system = FileReader.read_system_from(r'C:\Users\Pat\Documents\My Dropbox\Simulation Hecuba Group Share\ObjectOrientedOrbits\euler_problem.txt')
get_bodies_with_acceleration(system.bodies)


if __name__ == '__main__':
    print 'Running tests'
    unittest.main()

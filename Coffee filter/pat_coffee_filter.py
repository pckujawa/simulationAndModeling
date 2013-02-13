from __future__ import division
import math
import numpy as np
import pylab as pl
import collections

g = 9.8  # m/s^2

class CoffeeFilterFall(object):
    def __init__(self, initial_position, initial_velocity, terminal_velocity, method='linear'):
        '''method = 'linear' or 'quadratic' '''
        self.method = method
        self.terminal_velocity = terminal_velocity
        self.initial_position = initial_position
        self.initial_velocity = initial_velocity
        self.accelerations_at = collections.defaultdict(list)

    def step(self, t, state):
        v_prev = state[1]
        x_prev = state[0]
        v_value = v_prev / self.terminal_velocity
        if self.method != 'linear':
            v_value = v_value**2
        accel = g * (v_value - 1)
        # NOTE: don't do `self.accelerations.append(accel)` here because methods like Runge-Kutta eval more than once per step
        self.accelerations_at[t].append(accel)
        return np.array([v_prev, accel])

##    def run_once(self, integration_fn, elapsed_time, dt, current_state, predictor_corrector_prev_state=None):
##        '''Returns the next state'''
##        if predictor_corrector_prev_state is not None:
##            next_state = integration_fn(elapsed_time, predictor_corrector_prev_state, current_state, self.step, dt)
##        else:
##            next_state = integration_fn(elapsed_time, current_state, self.step, dt)
##        return next_state

    def run(self, integration_fn, dt, t_end=1, predictor_corrector_used=False):
        # Let us call the state of the system 'state'
        states = [np.array([self.initial_position, self.initial_velocity])]

        elapsed_time = 0
        while elapsed_time < t_end:
            if predictor_corrector_used:
                if len(states) == 1:
                    prev_state = None
                else:
                    prev_state = states[-2]
                next_state = integration_fn(elapsed_time, prev_state, states[-1], self.step, dt)
            else:
                next_state = integration_fn(elapsed_time, states[-1], self.step, dt)
            states.append(next_state)
            elapsed_time += dt

        # The trouble with the above, is that what you really want is all the positions
        # Turn the list to an array:
        states = np.array(states)

        # Now array slices get all positions and velocities:
        self.positions = states[:,0]
        self.velocities = states[:,1]
        self.times = pl.frange(0, elapsed_time, dt)

        # Ugly stuff to get the accelerations
        #  (certain methods, e.g. RK, take many steps at the same time value)
        # Ok, so actually there are 57 unique time points for acceleration. Of them, 27 are the same as the `times` and the rest are points about halfway in-between time points. It's a little too complicated to figure out an appropriate accel value for each time point (since they are floats and don't index well), so eff it for now.
##        print 'accelerations_at', self.accelerations_at
##        assert len(self.accelerations_at.keys()) == len(self.times), "Uh-oh, found different time points for acceleration."
##        self.accelerations = [sum(x)/len(x) for x in self.accelerations_at.itervalues()]

        # store vars in the object for easier inspection
        self.elapsed_time = elapsed_time
        self.states = states
##        l = locals().copy()
##        del l['self']
##        for key,value in l.iteritems():
##            setattr(self, key, value)
        return self  # allow for chaining (or being chained from ctor)


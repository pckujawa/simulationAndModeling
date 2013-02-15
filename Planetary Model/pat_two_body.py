from __future__ import division
import numpy as np
import pylab as pl

class MyTwoBody(object):
    def __init__(self, initial_state, step_fn):
        self.inititial_state = initial_state
        self.step_fn = step_fn

    def step(self, t, state):
        return self.step_fn(state, t)  # backwards param order

    def run(self, integration_fn, dt, t_end=1, predictor_corrector_used=False):
        # Let us call the state of the system 'state'
        states = [np.array(self.inititial_state)]

        elapsed_time = 0
        while elapsed_time < t_end:
            next_state = integration_fn(elapsed_time, states[-1], self.step, dt)
            states.append(next_state)
            elapsed_time += dt

        # The trouble with the above, is that what you really want is all the positions
        # Turn the list to an array:
        states = np.array(states)

        # Now array slices get all positions and velocities:
        # x1,y1,vx1,vy1,x2,y2,vx2,vy2 = state
        self.x1s  = states[:,0]
        self.y1s  = states[:,1]
        self.vx1s = states[:,2]
        self.vy1s = states[:,3]
        self.x2s  = states[:,4]
        self.y2s  = states[:,5]
        self.vx2s = states[:,6]
        self.vy2s = states[:,7]
        self.times = pl.frange(0, elapsed_time, dt)

        # store vars in the object for easier inspection
        self.elapsed_time = elapsed_time
        self.states = states
##        l = locals().copy()
##        del l['self']
##        for key,value in l.iteritems():
##            setattr(self, key, value)
        return self  # allow for chaining (or being chained from ctor)


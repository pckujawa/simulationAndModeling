#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
from scipy.integrate import odeint  # for integrate.odeint
import numpy as np
import pylab as pl
import math
import pat_ode_solvers as myode
from pat_two_body import MyTwoBody

# Notice that parameters are global.
GM = 4 * math.pi**2
ratio_m1_M    = 0.001
ratio_m2_M    = 0.04

def two_body(state, t):  #OUCH! The signature is reversed for odeint!
    # s_1 = (x1, y1, vx1, vy1), s_2 = (x2, y2, vx2, vy1)
    x1,y1,vx1,vy1,x2,y2,vx2,vy2 = state
    r1_sqrt_cube  = (x1**2 + y1**2)**(1.5)
    r2_sqrt_cube  = (x2**2 + y2**2)**(1.5)
    r21_sqrt_cube = ((x2-x1)**2 + (y2-y1)**2)**(1.5)
    a_1x = -GM * (x1/r1_sqrt_cube - ratio_m2_M*(x2-x1)/r21_sqrt_cube)
    a_1y = -GM * (y1/r1_sqrt_cube - ratio_m2_M*(y2-y1)/r21_sqrt_cube)
    a_2x = -GM * (x2/r2_sqrt_cube + ratio_m1_M*(x2-x1)/r21_sqrt_cube)
    a_2y = -GM * (y2/r2_sqrt_cube + ratio_m1_M*(y2-y1)/r21_sqrt_cube)
    return np.array([vx1, vy1, a_1x, a_1y, vx2, vy2, a_2x, a_2y])


class TwoBodyResult(object):
    def __init__(self, states):
        self.x1s  = states[:,0]
        self.y1s  = states[:,1]
        self.vx1s = states[:,2]
        self.vy1s = states[:,3]
        self.x2s  = states[:,4]
        self.y2s  = states[:,5]
        self.vx2s = states[:,6]
        self.vy2s = states[:,7]

## Problem setup
t_end = 100  # years
step_cnt = 1000
times = np.linspace(0.0, t_end, step_cnt)
yinit = np.array(
    [2.52, 0, 0, math.sqrt(GM/2.52), 5.24, 0, 0, math.sqrt(GM/5.24)])  # initial values

## My custom ODE
dt = t_end / step_cnt
my_ode = MyTwoBody(yinit, two_body)\
    .run(myode.runge_kutta, dt, t_end)
my_ode_results = TwoBodyResult(my_ode.states)

## Python ODE
states = odeint(two_body, yinit, times)
scipy_ode_results = TwoBodyResult(states)

# Plot the results

def plot_results(problem):
    pl.plot(problem.x1s, problem.y1s, 'bo-', label='m1')
    pl.plot(problem.x2s, problem.y2s, 'ro-', label='m2')
    pl.legend(loc='best')

pl.subplot(2,1,1)
pl.title('My ODE, RK, $dt=%.3f$' % dt)
plot_results(my_ode_results)
pl.subplot(2,1,2)
pl.title('SciPy ODE, $dt=%.3f$' % dt)
plot_results(scipy_ode_results)
pl.show()
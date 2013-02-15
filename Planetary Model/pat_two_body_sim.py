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

t_end = 100  # years
times = np.linspace(0.0, t_end, 1000)
yinit = np.array(
    [2.52, 0, 0, math.sqrt(GM/2.52), 5.24, 0, 0, math.sqrt(GM/5.24)])  # initial values

## My custom ODE
dt = 1
problem = MyTwoBody(yinit, two_body)\
    .run(myode.runge_kutta, dt, t_end)

pl.plot(problem.x1s, problem.y1s, 'bo-', label='m1')
pl.plot(problem.x2s, problem.y2s, 'ro-', label='m2')
pl.show()

## Python ODE
y = odeint(two_body, yinit, times)

# Plot the results
clf()
plot(y[:,0], y[:,1]) # y[:,0] is the first column of y, the postitions
xlabel('Times (s)')
ylabel('Positions (m)')
title('Two Body Motion')
grid()
show()

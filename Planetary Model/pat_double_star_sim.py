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
from matplotlib.animation import FuncAnimation  # v1.1+
import math
import pat_ode_solvers as myode


class Struct(object):
    '''Use as an expandable object, e.g. `x = Struct(); x.prop1 = "property"`'''
    pass

# Notice that parameters are global.
GM = 4 * math.pi**2

class Star():
    x = 0
    y = 0

star1, star2 = Star(), Star()
star1.x, star1.y = 0, 0
star2.x, star2.y = 2, 0

def double_star(state, t):
    '''Modify the program written for two planets to simulate the double star system, with the first star located at the origin and the second star of equal mass located at (2,0). Place the planet at (1.1,1) and systematically vary the x and y components of the velocity to obtain different types of orbits.'''
    x,y,vx,vy = state
    r1 = np.array([x - star1.x, y - star1.y])
    r2 = np.array([x - star2.x, y - star2.y])
    lenr1 = math.sqrt(r1.dot(r1))
    lenr2 = math.sqrt(r2.dot(r2))
    ax, ay = -GM * (r1/lenr1**3 + r2/lenr2**3)
    return np.array([vx, vy, ax, ay])


class DoubleStarResult(object):
    def __init__(self, states):
        self.xs  = states[:,0]
        self.ys  = states[:,1]
        self.vxs = states[:,2]
        self.vys = states[:,3]

## Problem setup
t_end = 100.0  # years
step_cnt = 1000
dt = t_end / step_cnt
times = np.linspace(0.0, t_end, step_cnt)

# Double star
x0, y0 = 1.1, 1
vx0, vy0 = math.sqrt(GM), math.sqrt(GM)
double_star_state_0 = np.array([x0, y0, vx0, vy0])
double_star_analysis = Struct()  # to hold all results

## Python ODE
states = odeint(double_star, double_star_state_0, times)
double_star_analysis.scipy_results = DoubleStarResult(states)


def do_plots():
    results = double_star_analysis.scipy_results
##    fig = pl.figure(figsize=(14,8))
##    fig.subplots_adjust(bottom=0.03, left=0.04, top = 0.95, right=0.95)
##    pl.title('SciPy ODE, $dt={:.3f}$, a/rtol={}/{}'.format(dt, atol, rtol))
    pl.title('SciPy ODE, $dt={:.3f}$'.format(dt))
    pl.ylabel('AU')
    pl.xlabel('AU')
    pl.scatter([star1.x], [star1.y], color='red')
##    pl.annotate('star 1',
##         xy=(star1.x, star1.y), xycoords='data',
##         xytext=(-60, 40), textcoords='offset points', fontsize=16,
##         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.5"))

    pl.scatter([star2.x], [star2.y], color='blue')
##    pl.annotate('star 2',
##         xy=(star2.x, star2.y), xycoords='data',
##         xytext=(-60, 40), textcoords='offset points', fontsize=16,
##         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.5"))

    pl.plot(results.xs, results.ys, 'o-', color='black', markersize=1, linewidth=0.1)

    info = '''Initial Conditions
$x_0={}$
$y_0={}$
$v_{{x0}}={:.3f}$
$v_{{y0}}={:.3f}$'''.format(x0, y0, vx0, vy0)
    pl.annotate(info, xy=(0.01, 0.01), xycoords='axes fraction', fontsize=12)

##    pl.legend(loc='best')
##    ax.yaxis.set_major_locator(pl.MaxNLocator(nbins=4))
##    pl.savefig('pat_orbits_energies_momentums_atol={}_rtol={}.png'.format(
##        atol, rtol))
    pl.show()

# Animate orbit
# Code courtesy of George Lesica
results = double_star_analysis.scipy_results
fig = pl.figure(figsize=(8, 8))
ax = pl.axes(xlim=(-6, 6), ylim=(-6, 6))
pl.title('Orbit Simulation', fontsize=16)
pl.xlabel('X Position')
pl.ylabel('Y Position')
p1, = ax.plot([], [], marker='o')
p2, = ax.plot([], [], marker='o')
p3, = ax.plot([], [], marker='o')

def animate(i):
    p1.set_data([star1.x], [star1.y])
    p2.set_data([star2.x], [star2.y])
    p3.set_data([results.xs[i]], [results.ys[i]])
    return p1, p2, p3

anim = FuncAnimation(fig, animate, frames=1000, interval=1, blit=True)
anim.save('pat_double_star_animation.mp4', fps=30)
pl.show()

##do_plots()


## Animation (key steps from Jesse's code)
##pl.ion()
##plot_handle = pl.plot(results.xs, results.ys, 'o-', color='black', markersize=1, linewidth=0.1)[0]
##for...
##    plot_handle.set_data()
##pl.ioff()
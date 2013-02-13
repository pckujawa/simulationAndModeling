#-------------------------------------------------------------------------------
# Purpose:
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python

from __future__ import division
import math
import numpy as np
import pylab as pl
from pprint import pprint, pformat
import pat_ode_solvers as ode
import pat_coffee_filter as cf

empirical_ts_ys = np.array([[.2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,.4280,
              .4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,.6498,.6744,
              .6990,.7236,.7482,.7728,.7974,.8220,.8466],
              [.4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,.3609,
               .3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,.2497,.2337,
               .2175,.2008,.1846,.1696,.1566,.1393,.1263]])

emp_positions = empirical_ts_ys[1,:]  # meters
emp_times     = empirical_ts_ys[0,:]  # seconds
emp_times -= emp_times[0]  # shift to zero on the first measurement


# Take a look at about the area in which velocity seems constant (terminal)
v_term_points = emp_times > 0.3
ys_term = emp_positions[v_term_points]
times_term = emp_times[v_term_points]
v_term = (ys_term[-1] - ys_term[0]) / (times_term[-1] - times_term[0])  # negative
v_term_times     = np.array([emp_times[v_term_points][0], emp_times[-1]])
v_term_endpoints = ys_term[0] * np.array([1, 1 + v_term])

## Finite differencing of emp data to guess velocity and accel
dts = emp_times[1:] - emp_times[:-1]  # 26 items
dts = dts[:-1]  # ditch last dt (seems more appropriate than first dt)
fd_vs = (emp_positions[2:] - emp_positions[:-2]) / (2*dts)
fd_as =(emp_positions[2:] + emp_positions[:-2] - 2*emp_positions[1:-1]) / (dts**2)


## Simulation
y_0, v_0 = emp_positions[0], 0
t_end = emp_times[-1]  # relies on times being shifted
dt = dts.mean()
linear_sim = cf.CoffeeFilterFall(y_0, v_0, v_term, method='linear')\
    .run(ode.runge_kutta, dt, t_end = t_end)
quadratic_sim = cf.CoffeeFilterFall(y_0, v_0, v_term, method='quadratic')\
    .run(ode.runge_kutta, dt, t_end = t_end)


def make_plot_better():
    # Move axes to be on the plot
    ax = pl.gca()  # get current axes
    ax.xaxis.set_major_locator(pl.MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(pl.MaxNLocator(nbins=4))

## Plotting
sims = [(linear_sim, 'r'), (quadratic_sim, 'k')]
pl.figure(figsize=(12, 9))
pl.subplot(4,1,1)
pl.title('Coffee Filter Analyses')
pl.ylabel('Position (m)')
pl.plot(v_term_times, v_term_endpoints, label='$v_{term}$ guess')
pl.scatter(emp_times, emp_positions, label='empirical')
for sim,c in sims:
    pl.scatter(sim.times, sim.positions, c=c, label=sim.method)
pl.xlim(-0.2, 0.7)
make_plot_better()
pl.legend(loc='lower left', frameon=False)  #TODO make smaller/alpha

pl.subplot(4,1,2)
pl.ylabel('Velocity (m/s)')
fd_times = emp_times[0] + np.cumsum(dts)
pl.scatter(fd_times, fd_vs, label='finite diff')
for sim,c in sims:
    pl.scatter(sim.times, sim.velocities, c=c, label=sim.method)
##    pl.plot(sim.times, sim.velocities)
pl.xlim(-0.2, 0.7)
pl.plot(fd_times, [v_term]*len(fd_times), label=('$v_{term}=%.2f$' % v_term))
make_plot_better()
pl.legend(loc='lower left', frameon=False)  #TODO make smaller/alpha

pl.subplot(4,1,3)
pl.ylabel(r'Acceleration (m/$s^2$)')
pl.scatter(fd_times, fd_as, label='finite diff')
pl.xlim(0, emp_times[-1])
pl.xlabel('Time (s)')
pl.subplots_adjust(hspace=0.3)
# SIM accel is too hard
##for sim in sims:
##    pl.scatter(sim.times, sim.accelerations, 0.5, label=sim.method)
make_plot_better()
pl.xlim(-0.2, 0.7)

pl.subplot(4,1,4)
pl.ylabel('Acceleration (m/$s^2$)')
pl.xlabel('Velocity (m/s)')
pl.scatter(fd_vs, fd_as, label='fd')
pl.plot([fd_vs.min(), fd_vs.max()], [0, 0], 'b-')  # line at y=0
# SIM accel is too hard
##for sim in sims:
##    pl.scatter(sim.velocities, sim.accelerations, 0.5, label=sim.method)
make_plot_better()
pl.savefig("pat_coffee_filter_sim.png", dpi=72)
pl.show()

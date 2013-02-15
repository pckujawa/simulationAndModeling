#-------------------------------------------------------------------------------
# Purpose:
# Author:      Pat
#-------------------------------------------------------------------------------
#!/usr/bin/env python

from __future__ import division

def euler(t, x, f, dt):
    return x + f(t, x) * dt

def euler_richardson_midpoint(t, x, f, dt):
    halfdt = 0.5 * dt
    xmid = x + halfdt * f(t, x)
    xnext = x + f(t + halfdt, xmid)*dt
    return xnext

def runge_kutta(t, x, f, dt):
    k1 = f(t         , x         )*dt
    k2 = f(t + 0.5*dt, x + 0.5*k1)*dt
    k3 = f(t + 0.5*dt, x + 0.5*k2)*dt
    k4 = f(t + dt    , x + k3    )*dt
    xnext = x + 1/6.0*(k1 + 2*k2 + 2*k3 + k4)
    return xnext

def predictor_corrector(t, x_prev, x, f, dt):
    # Take first step with RK to init
    if x_prev is None:
        x_prev = runge_kutta(t, x, f, dt)
    x_p = x_prev + 2*f(t, x) * dt  # _predictor
    xnext = x + 0.5 * (f(t, x_p) + f(t, x)) * dt
    return xnext

if __name__ == '__main__':
    print 'no tests, only use as a library'

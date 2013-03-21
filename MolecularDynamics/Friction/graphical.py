#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction graphics
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import pylab as pl
##import math
##import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+
import os

working_directory = os.getcwd()

# Circle code courtesy of Kevin Joyce
def get_nice_circle(x, y, radius, color="lightsteelblue", facecolor="green", alpha=0.6, ax=None):
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


def animate_with_live_integration(
    containers, integrator, dt, xlim, ylim, figsize, particle_radius,
    frame_show_modulus=1, num_frames_to_bootstrap=0, info_for_naming='', save_animation=False, show_animation=True, run_until_func=None):
    """
    """
    global num_forward_frames
    init_container = containers[0]
    plot_title = 'Molec Dyn Friction' + info_for_naming

    also_run_backwards = False
    num_forward_frames = len(containers)

    fig = pl.figure(figsize=figsize)
    ax = pl.axes(xlim=xlim, ylim=ylim)
    ax.set_aspect('equal')  # NOTE: This can make the plot look squished if figsize isn't wide enough
    pl.title(plot_title, fontsize=16)
    pl.xlabel('X Position')
    pl.ylabel('Y Position')

    circles = []
    for i in xrange(init_container.num_particles):
        e = get_nice_circle(0, 0, 0.5*particle_radius)
        circles.append(ax.add_patch(e))
    time_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    pulling_force_text = ax.text(0.5, 0.90, '', transform=ax.transAxes)
    damp_force_text = ax.text(0.5, 0.80, '', transform=ax.transAxes)

    # Bootstrap some frames
    if run_until_func is not None:
        while run_until_func(containers[-1]):
            containers.append(integrator.step(containers[-1], dt))
    num_forward_frames = len(containers) // frame_show_modulus  # NOTE: int division
    while num_forward_frames < num_frames_to_bootstrap:
        num_forward_frames += 1
        for i in xrange(frame_show_modulus):
            next_container = integrator.step(containers[-1], dt)
            containers.append(next_container)

    # init fn seems to prevent 'ghosting' of first-plotted data
    def init():
        """initialize animation"""
        time_text.set_text('')
        pulling_force_text.set_text('')
        damp_force_text.set_text('')
        for c in circles:
            c.center = (-1, -1)  # hide off-screen
        return circles + [time_text, pulling_force_text, damp_force_text]


    def next_frame(ix_frame):
        global num_forward_frames
        if ix_frame % 100 == 1:
            print 'frame', ix_frame
        # Do only the integration necessary to get to the requested frame
        container_ix = ix_frame*frame_show_modulus
        while num_forward_frames <= ix_frame:
            num_forward_frames += 1
            for _ in xrange(1 + frame_show_modulus):  # always run at least once
                next_container = integrator.step(containers[-1], dt)
                containers.append(next_container)
        c = containers[container_ix]
        posns = c.positions
        try:
            time_text.set_text('time = %.1f' % c.time)
            pulling_force_text.set_text('Fp = ' + str(c.pull_accelerations))
            damp_force_text.set_text('Fd = ' + str(c.damp_accelerations))
        except AttributeError:
            pass
        facecolor = 'green'
        # TODO paint floor black, sled green, etc
        for i,circle in zip(xrange(init_container.num_particles), circles):
            circle.center = (posns[i][0], posns[i][1])  # x and y
            circle.set_facecolor(facecolor)
        return circles + [time_text, pulling_force_text, damp_force_text]

    # Amount of framedata to keep around for saving movies.
    save_count = 5000  # should be enough for animation to be able to save entire movie, if desired, without having to re-run
    anim = FuncAnimation(fig, next_frame, interval=dt, blit=True, init_func=init, save_count=save_count)
    if show_animation:
        pl.show()

    end_container = containers[-1]

##    # Now run... backwards!
##    if also_run_backwards:
##        for i in xrange(num_forward_frames):
##            next_container = integrator.step(containers[-1], -dt)
##            containers.append(next_container)

    num_total_frames = num_forward_frames
    if also_run_backwards:
        num_total_frames += num_forward_frames

    if save_animation:
        if not show_animation:
            # If we haven't run/showed yet, we need to have the animation know how many frames to run for
            anim = FuncAnimation(fig, next_frame, interval=dt, blit=True, frames=num_total_frames)
        print 'beginning save of animation with frame count:',   num_total_frames
        try:
            anim.save('{}/anim/{}.avi'.format(working_directory, info_for_naming), fps=30, dpi=72)
            pl.close()
        except:  # Tk error, although it still seems to bubble up
            pass


def plot_potential_energy(containers, dt, num_total_frames):
    init_container = containers[0]
    pl.clf()
    # Plot PE --------------------------
    times = pl.frange(dt, num_total_frames*dt, dt)  # skip zeroth time because we have no value for it
    pes = [c.potential_energy for c in containers[1:]]  # skip first container because it has no PE
    plotted_pes = np.array(pes[:len(times)])
    plotted_pes /= init_container.num_particles  # PE per particle
    pl.plot(times, plotted_pes, 'o-', color='black', markersize=1, linewidth=0.1)
    pl.ylabel('Potential energy per particle')
    pl.xlabel('Time')
    pl.title('PE/particle for {} frames'.format(num_forward_frames))
    pl.savefig('{}/plots/pe {}.png'.format(working_directory, info_for_naming))
    pl.show()


def plot_pulling_force(times, pulling_forces, info_for_naming=''):
    fig = pl.figure()
    pulling_forces = np.array(pulling_forces)
    fx_max_ix = np.argmax(pulling_forces, axis=0)[0]
    fx_max = pulling_forces[fx_max_ix][0]
    time_of_fx_max = times[fx_max_ix]
    defaults = {'markersize':1, 'linewidth':0.1}
    for ix_dim, name in zip(range(len(pulling_forces[0])), ['Fx']):#, 'Fy', 'Fz']):
        pl.plot(times, pulling_forces[:, ix_dim], 'o-', label=name, **defaults)
    pl.ylabel('Fp')
    pl.xlabel('Time')
    pl.title('Pulling force vs time')
    # Max Fx:
    pl.annotate('$F_x={:.1f}$'.format(fx_max), xy=(time_of_fx_max, fx_max),
        xytext=(0.2, 0.9), textcoords='axes fraction', fontsize=14,
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.1"))
    # Sim info:
    pl.ylim(ymin=-2)  # make room for it
    pl.annotate(info_for_naming, xy=(0.01, 0.03), xycoords='axes fraction', fontsize=12)
    pl.savefig('{}/plots/Fp {}.png'.format(working_directory, info_for_naming))
    pl.show()


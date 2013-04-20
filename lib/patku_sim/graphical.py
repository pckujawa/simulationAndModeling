#-------------------------------------------------------------------------------
# Name:        molecular dynamics friction graphics
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np
import pylab as pl
import math
##import unittest
##import itertools
##from pprint import pprint, pformat
##from collections import defaultdict, namedtuple
from matplotlib.animation import FuncAnimation  # v1.1+
import os
import textwrap
import time
from libs import make_dirs_if_necessary

working_directory = os.getcwd()

# Circle code courtesy of Kevin Joyce
def get_nice_circle(x, y, radius, color="black", facecolor="tan", alpha=0.6, ax=None):
    """ add a circle to ax or current axes
    """
    e = pl.Circle([x, y], radius)
    if ax is None:
        ax = pl.gca()
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_edgecolor( color )
    e.set_linewidth(0)
    e.set_facecolor( facecolor )  # "none" not None
    e.set_alpha( alpha )
    return e


def animate_with_live_integration(
    containers, integrator, dt, xlim, ylim, figsize, particle_radius,
    frame_show_modulus=1, num_frames_to_bootstrap=0, info_for_naming='', save_animation=False, show_animation=True, run_until_func=None, anim_save_kwargs=None, anchor_ixs=None, sim_wide_params=None):
    """
    """
    global num_forward_frames

    anim_save_kwargs = anim_save_kwargs or {'fps': 30}
    init_container = containers[0]
    plot_title = info_for_naming

    also_run_backwards = False
    num_forward_frames = len(containers)

    for i in [0, 1]:
        assert type(figsize[i]) is int, "Can't save an animation with float figsize"
    fig = pl.figure(figsize=figsize)
    ax = pl.axes(xlim=xlim, ylim=ylim)
    ax.set_aspect('equal')  # NOTE: This can make the plot look squished if figsize isn't wide enough
    title_width = int(figsize[0] * 10)  # dynamic way to figure out how much text will fit in the title
    title_wrapped = '\n'.join(textwrap.wrap(plot_title, width=title_width))
    pl.title(title_wrapped, fontsize=12)
##    pl.xlabel('X Position')
##    pl.ylabel('Y Position')

    # Bootstrap some frames using user-suppliend func and/or count
    if run_until_func is not None:
        while run_until_func(containers[-1]):
            containers.append(integrator.step(containers[-1], dt))
    num_forward_frames = len(containers) // frame_show_modulus  # NOTE: int division
    while num_forward_frames < num_frames_to_bootstrap:
        num_forward_frames += 1
        for i in xrange(frame_show_modulus):
            next_container = integrator.step(containers[-1], dt)
            containers.append(next_container)

    ## Set up plotting
    circles = []
    for i in xrange(init_container.num_particles):
        e = get_nice_circle(0, 0, particle_radius)
        circles.append(ax.add_patch(e))

    # NOTE: We can animate an annotations, but not using textcoords='figure points'
    multi_text = ax.annotate('stats', (0, 0),
            xytext=(5, 5), textcoords='axes points',
            ha='left', va='bottom')
    time_text = ax.text(0, -2, 'time', transform=ax.transAxes)
    flux_text = ax.text(0.5, -2, 'flux', transform=ax.transAxes)
    pulling_force_text = ax.text(0.5, 0.90, '', transform=ax.transAxes)
    damp_force_text = ax.text(0.5, 0.80, '', transform=ax.transAxes)
##    aperture_line = ax.plot(xlim, [sim_wide_params.y_funnel_bottom]*2,
##            linestyle='-', color='black', alpha=0.5,
##            linewidth=0.5)

    # init fn seems to prevent 'ghosting' of first-plotted data
    texts = [multi_text, time_text, flux_text, pulling_force_text, damp_force_text]
    def init():
        """initialize animation"""
        for t in texts:
            t.set_text('')
        for c in circles:
            c.center = (-100, -100)  # hide off-screen
        return circles + texts


    def next_frame(ix_frame):
        global num_forward_frames
##        if ix_frame % 100 == 1:
##            print 'frame', ix_frame
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
            stats = 'time={t:.1f}\ncount={n}'.format(
                    t = c.time, n = c.cumulative_grains_below_aperture)
            multi_text.set_text(stats)
##            time_text.set_text('time = %.1f' % c.time)
##            flux_text.set_text('count = {}'.format(c.cumulative_grains_below_aperture))
            pulling_force_text.set_text('Fp = ' + str(c.pull_accelerations))
            damp_force_text.set_text('Fd = ' + str(c.damp_accelerations))
        except AttributeError: pass
        for i, circle in zip(xrange(init_container.num_particles), circles):
            circle.center = (posns[i][0], posns[i][1])  # x and y
        try:
            for ix, force_mag in zip(sim_wide_params.anchor_ixs, c.anchor_accels):
                anchor_circle = circles[ix]
                anchor_circle.set_facecolor('black')
                # Consider using logarithm
                max_expected_force = 1000
                alpha = min(math.log(force_mag + 2, max_expected_force), 1)  # always > 0.1 so visible; cap at 100%
##                alpha = min(force_mag/100.0 + 0.01, 1)
                anchor_circle.set_alpha(alpha)
        except AttributeError: pass
        return circles + texts

    # Amount of framedata to keep around for saving movies.
    save_count = 5000  # should be enough for animation to be able to save entire movie, if desired, without having to re-run
    anim = FuncAnimation(fig, next_frame,
            interval=dt, blit=True, init_func=init, save_count=save_count)
    if show_animation:
        try:
            pl.show()
        except AttributeError as e:  # 'NoneType' object has no attribute 'tk'
            print e
    else:
        pl.close()

##    end_container = containers[-1]
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
        print 'beginning save of animation with frame count:',   num_total_frames, 'at', time.strftime("%a %b %d %H:%M:%S")
        prefix = '{wd}/dumps/{name}'.format(wd = working_directory,
                name = info_for_naming)
        make_dirs_if_necessary(prefix)
        anim.save('{pre}/anim len={len}.avi'.format(
                pre=prefix, len=num_total_frames), **anim_save_kwargs)
        print 'done saving at', time.strftime("%a %b %d %H:%M:%S")
        try:
            pl.close()
        except:
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
    pl.savefig('{}/plots/svg/pe {}.svg'.format(working_directory, info_for_naming))
    pl.savefig('{}/plots/pe {}.png'.format(working_directory, info_for_naming))
    pl.show()


def get_pulling_force_max_and_meta(times, pulling_forces, normal_force=0):
    """
    """
    drop_magnitude = 0.01 * (normal_force + 90)  # roughly empirical. Want at most 0.9 at W=-20 but > 1 at W=40, based off plots
    fx_max = 0
    for ix, (Fp, time) in enumerate( zip(pulling_forces[:, 0], times) ):
        if Fp < fx_max - drop_magnitude:  # if we get a significant drop, call that max
            break
        if Fp > fx_max:
            fx_max = Fp
            time_of_fx_max = time
            fx_max_ix = ix
##    fx_max_ix = np.argmax(pulling_forces, axis=0)[0]
##    fx_max = pulling_forces[fx_max_ix][0]
##    time_of_fx_max = times[fx_max_ix]
    print 'Fx_max of {} at t={}'.format(fx_max, time_of_fx_max)
    return fx_max, time_of_fx_max, fx_max_ix

def plot_pulling_force(times, pulling_forces, normal_force, info_for_naming='', show=True):
    pulling_forces = np.array(pulling_forces)  # may be redundant
    fx_max, time_of_fx_max, fx_max_ix = get_pulling_force_max_and_meta(times, pulling_forces, normal_force=normal_force)
    defaults = {'markersize':1, 'linewidth':0.1}
    fig = pl.figure()
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
    pl.savefig('{}/plots/svg/Fp {}.svg'.format(working_directory, info_for_naming))
    pl.savefig('{}/plots/Fp {}.png'.format(working_directory, info_for_naming))
    if show:
        pl.show()
    else:
        pl.close()


def plot_all_pulling_forces(lstats, filename='all', show=True):
    defaults = {'markersize':1, 'linewidth':1, 'alpha': 0.5}
    colors = {1: 'purple', 9: 'black', 13: 'red', 17: 'blue'}  # by sled size
    fig = pl.figure(figsize=(14,10))
    for stats in lstats:
        label = '{} {:3d} {:.1f}'.format(
            stats.sled_size, stats.W, stats.Fp_max)
        pulling_forces = stats.pulling_forces
        times = stats.times
        for ix_dim, dim in zip(range(len(pulling_forces[0])), ['Fx']):#, 'Fy', 'Fz']):
            pl.plot(times, pulling_forces[:, ix_dim], '-',
                label=label, color=colors[stats.sled_size],
                **defaults)
    pl.ylabel('Fp')
    pl.xlabel('Time')
    pl.title('Pulling force vs time, varying sled size and normal force')
##    pl.annotate('{}'.format(colors), xy=(0.1, 0.9),
##        xycoords='axes fraction', fontsize=12)
##    pl.xlim(xmin=-10)  # make room for legend
##    ax = pl.gca()
##    handles, labels = ax.get_legend_handles_labels()
##    pl.legend(handles, labels)
    pl.legend(ncol=3, title='Sled size, normal force, max Fp',
        fontsize=10, loc='upper left', frameon=False)
##    pl.savefig('{}/plots/svg/Fp {}.svg'.format(working_directory, filename))  # SVG is ~30MB and takes a long time
    pl.savefig('{}/plots/Fp {}.png'.format(working_directory, filename))
    if show:
        pl.show()
    else:
        pl.close()
##plot_all_pulling_forces(data_table)


def plot_friction_slope(data_frame, filename='friction slope', show=True):
    defaults = {'markersize':5, 'linewidth':1, 'alpha': 0.5}
    colors = {1: 'purple', 9: 'black', 13: 'red', 17: 'blue'}  # by sled size
    fig = pl.figure(figsize=(8,6))
    for sled_size in np.unique(data_frame['sled_size']):
        aggregate = data_frame[data_frame['sled_size'] == sled_size]
        x = aggregate['W']
        y = aggregate['Fp_max']
        slope, yint = np.polyfit(x, y, deg=1)
        area = (sled_size + 1) / 2  # just the floor of the sled
        c = yint / area  # constant multiplier
        label = 'Sled of {:.0f}: $f_s = {:.3f}W + {:.3f}A$'.format(sled_size, slope, c)
        color = colors[sled_size]
        pl.plot(x, y, 'o', color=color, **defaults)
        xfit = np.array([min(x), max(x)])
        yfit = slope * xfit + yint
        pl.plot(xfit, yfit,
            '-', label=label, color=color, **defaults)
    xmin, xmax = pl.xlim()
    pl.xlim((xmin-1, xmax+1))
    pl.ylabel('$max(Fp)$')
    pl.xlabel('W')
    pl.title('Friction force (max pulling force vs normal force)')
    ax = pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    pl.legend(handles[::-1], labels[::-1], # reverse order
        ncol=1, fontsize=14, loc='upper left', frameon=False)
    pl.savefig('{}/plots/svg/Fp {}.svg'.format(working_directory, filename))
    pl.savefig('{}/plots/Fp {}.png'.format(working_directory, filename))
    if show:
        pl.show()
    else:
        pl.close()
##plot_friction_slope(data_frame)

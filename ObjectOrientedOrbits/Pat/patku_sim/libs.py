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
##from matplotlib.animation import FuncAnimation  # v1.1+


class FileReader(object):
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
        pass
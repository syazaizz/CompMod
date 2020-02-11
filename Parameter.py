"""
Project Solar System: Parameter, a class to describe parameters for 3D planets

Version: 01/2020
"""

import numpy as np

class Parameter3D(object):
    """
    Class to describe parameters for 3D planets.

    Properties:
    dt (float) - timestep used in the simulation
    numstep (int) - number of steps for the simulation

    Methods:
    * formatted output
    """

    def __init__(self, dt, numstep):
        """
        Initialise parameters for Particle3D instance

        :param a: curvature as float
        :param D_e: depth as float
        :param r_e: position as float
        """
        self.dt = dt
        self.numstep = numstep


    def __str__(self):
        """
        Define output format.
        For particle parameter this will print as
        "a D_e r_e"
        """
        return "dt = " + str(self.dt) +" "+ "numstep = " + str(self.numstep)


    #read content from filename
    # call Particle3D __init__ method
    @staticmethod
    def from_file(file_handle):
        line = file_handle.readline()
        data = line.split(",")
        dt = np.array([float(i) for i in data[0:1]]) #data[0]
        numstep = np.array([float(i) for i in data[1:2]])
        return Parameter3D(dt, numstep)

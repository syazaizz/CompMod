"""
Project Solar System: ParticleManyBody main module.
THe program uses Velocity Verlet time integration mehod to simulate the main bodies of our solar systems,
interacting through Newtonian gravity.

Version: 01/2020
"""

import sys
import math
import numpy as np
from Particle import Particle
from scipy.signal import argrelextrema as ext

# Begin main code

def count(particle_file):
    count = 0
    with particle_file as f:
        for line in f:
            count += 1
        return count

def main():

    # Read name of output file from command line
    if len(sys.argv)!=4:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <parameter_file>"+ " <particle_file>"+ "output_trajectory_file")
        quit()
    else:
        parameter_file = sys.argv[1]
        particle_file_name = sys.argv[2]
        outfile_trajectory = sys.argv[3]


    #parameter of the simulation
    #filename1.txt contain data for parameter
    param_file = open(parameter_file, "r")
    line = param_file.readline()
    data = line.split(",")
    dt = np.array([float(i) for i in data[0:1]]) #data[0]
    numstep = np.array([float(i) for i in data[1:2]])
    print("dt =", dt, "and numstep =", numstep)


    #set up particles initial condition #list n loops to read p1, p2 ...
    particle_file = open(particle_file_name, "r")
    count_1 = count(particle_file)
    particle_file = open(particle_file_name, "r")
    particle=[]
    for i in range(0,count_1):
         particle.append(Particle.from_file(particle_file))

    for pp in particle:
        print(pp)

    """
    particle = 1
    while particle != 0:v_o = particle.velocity
        particle = Particle.from_file(particle_file)
        print(particle)
    """


    #COM correction
    M = 0                     #whats diff put M=0 or pp.mass
    P = 0
    for pp in particle:
        print(pp.velocity)
        M += pp.mass
        P += pp.velocity*pp.mass
    v_com = P/M

    #define new particle velocity
    for pp in particle:
        pp.velocity -= v_com
        print(pp.velocity)


    # Open output file
    trajectory_file = open(outfile_trajectory, "w")

    #write to trajectory_file
    #trajectory_file.write(p1)    #loop many part


    # Write out initial conditions
    #Total Energy

main()

"""
Project: Solar System
ManyBody Module (Main)

Use of the Velocity Verlet time integration method to simulate main bodies of our Solar System,
interacting through Newtonian Gravity.

Version: JP/03/0220
"""

import sys
import math
import numpy as np
from Particle import Particle
from scipy.signal import argrelextrema as extrema


# METHOD: Count the number of particles presented in a text file.
def count(particle_file):
    count = 0
    with particle_file as f:
        for line in f:
            count += 1
    return count


# MAIN CODE
def main():

    # COMMAND LINE VERIFICATION: Read the arguments presented in the command line.
    # Ensure the correct number of arguments are presented.
    if len(sys.argv) != 4:
        print("Incorrect number of arguments; you need four.")
        print("Usage: " + sys.argv[0] + " <parameter file>" + " <particle file>" + " <output file>")
        quit()
    else:
        parameter_file = sys.argv[1]
        particle_file = sys.argv[2]
        outfile = sys.argv[3]
    
    # Open output file
    trajectory = open(outfile, "w")


    # SIMULATION PARAMETERS
    # Read in to program from file.
    parameters = open(parameter_file, "r")            # Opens parameter file
    line = parameters.readline()                      # Reads each line from file
    tokena = line.split(",")                          # Splits line into a list, using comma
    dt = tokena[0]                                    # Timestep read from first element
    numstep = tokena[1]                               # Range read from last element
    time = 0.0                                        # Initial time set to zero

    print("Timestep (days) = " + dt + "; Range (days) = " + numstep)


    # READ PARTICLES AND INITIAL CONDITIONS: Particle information and initial
    # conditions are loaded into lists that are operated on during the simulation.
    file_handle = open(particle_file, "r")
    quantity = count(file_handle)                     # This also closes the file; needs reopened
    file_handle = open(particle_file, "r")
    particles = []                                    # Opens empty list
    for i in range(0, quantity):
        particles.append(Particle.from_file(file_handle))    #list filled with Particle instances
    #print(quantity)


    # CORRECT INITIAL VELOCITIES BY CENTRE OF MASS MOMENTUM
    sysmomentum = np.array([0.0, 0.0, 0.0])                  # Assignment of momentum vector
    sysmass = 0.0                                            # Assignment of mass

    for object in particles:
        sysmomentum += ((object.velocity)*(object.mass))     # Calculates total momentum of system
        sysmass += object.mass                               # Calculates total mass in system

    com_correction = sysmomentum / sysmass                   # Calculates linear drift velocity

    for object in particles:
        #print(object.velocity)
        object.velocity -= com_correction                    # Initial velocities are adjusted
        #print(object.velocity)


    # ASSIGN VARIABLES FOR ENERGY
    Energy_Kinetic = 0.0
    Energy_Potential = 0.0
    

    # INITIAL FORCE CALCULATION (5.2)
    calc = np.zeros(3)
    forces = []
    #forces = np.zeros((quantity, 3))
    #print(forces)
    for i in particles:
        for j in particles:
            if i != j:
                calc += i.force(i,j)
        #print(calc)
        forces.append(calc)
        calc = np.zeros(3)
    print(forces)


    # TIME INTEGRATION LOOP
    time = 0
    for i in particles:

        # UPDATE PARTICLE POSITIONS
        a = 0
        b = 0
        #print(forces[1][2])
        for j in particles:
            j.leap_pos2nd(dt, float(forces[a]))
            a += 1
        
        # UPDATE FORCES
        calc = np.zeros(3)
        force_update = []
        for k in particles:
            for l in particles:
                if k != l:
                    calc += k.force(k,l)
            print(calc)
            force_update.append(calc)

        # UPDATE PARTICLE VELOCITY
        # Requires averaging the current and new forces
        a = 0
        for j in particles:
            j.leap_velocity(dt, 0.5*(forces[a]+force_update[a]))
            a += 1

        # REDEFINE FORCES WITH UPDATED VALUE
        for j in range(0, quantity):
            forces[j] = force_update[j]

        # INCREMENT TIME
        time += dt
        
        # WRITE DATA TO OUTPUT FILE
        trajectory.write(quantity)
        trajectory.write("\n")
        trajectory.write("Point = " + dt)
        trajectory.write("\n")
        for iterate in particles:
            trajectory.write(iterate.__str__())
            trajectory.write("\n")
        trajectory.close()


main()

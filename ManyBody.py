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
#from itertools import cycle




#count number of particle in particle.txt file
def count(particle_file):
    count = 0
    with particle_file as f:
        for line in f:
            count += 1
        return count

# Begin main code
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

    # Open output file
    trajectory_file = open(outfile_trajectory, "w")


    #parameter of the simulation
    #filename1.txt contain data for parameter
    param_file = open(parameter_file, "r")
    line = param_file.readline()
    data = line.split(",")
    dt = np.array([float(i) for i in data[0:1]]) #data[0]
    numstep = np.array([float(i) for i in data[1:2]])
    time = 0.0
    print("dt =", dt, "and numstep =", numstep)


    #set up particles initial condition #list n loops to read p1, p2 ...
    particle_file = open(particle_file_name, "r")
    no_particle = count(particle_file)
    particle_file = open(particle_file_name, "r")
    particle=[]
    for i in range(0,no_particle):                                    #force
         particle.append(Particle.from_file(particle_file))

    # for pp in particle:
    #     print(pp)


    #COM correction
    M = 0                               #whats diff put M=0 or pp.mass?
    P = 0
    for pp in particle:
        #print(pp.velocity)
        #M = 0                               #pos M and P ?
        #P = 0
        M += pp.mass
        P += pp.velocity*pp.mass
    v_com = P/M
    #print(v_com)

    KE = 0
    PE = 0
    force = [0, no part]
    pos_list = []
    sep_list = []                                                #might need to specify sep between sun, and moon earth?

    for pp in particle:
        #define new particle velocity
        pp.velocity -= v_com
        #print(pp.velocity)


        #check: https://www.mdanalysis.org/docs/_modules/MDAnalysis/coordinates/XYZ.html            #loops
        def VMD(no_particle, time, pp):
            return str(no_particle) + '\n' + str(time) + '\n' + str(pp)                             #check pp?? VMD format
            #return str("{0:d}\n".format(no_particle)) + str("{0}\n".format(time)) + str(pp)        #whats this?

        VMD = VMD(no_particle, time, pp)
        print(VMD)


        #write to trajectory_file
        trajectory_file.write(VMD)

        pos = np.linalg.norm(pp.position)
        pos_list.append(pos)
        #print(pos_list)

        sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))                            #how to exclude the 1st row?
        sep_list.append(sep_pos)


        # Write out initial energy and force
        KE += Particle.kin_energy(pp)


        for p in particle:

            if pp != p:

                PE += Particle.pot_energy(pp,p)/2

                force += Particle.force(pp,p)     #cxheck unit G

    # print(KE)
    # print(PE)
    # print(force)


    #initial total energy
    energy = KE + PE
    print(energy)                           #check: total energy - ignore -ve sign PE? KE + PE got -ve energy, KE-PE: +ve energy


    #initialise data lists
    time_list = [time]
    energy_list = [energy]
    # print(pos_list)
    # print(sep_list)

    force_new = 0

    for pp in particle:
        for p in particle:

            for i in range(int(numstep)):  #outermost

                # Update particle position
                pp.leap_pos2nd(dt, force)             #wrong
                p.leap_pos2nd(dt, force)

                # Update force
                if pp != p:
                    force_new += Particle.force(pp,p)   #fx
                #print(force_new)

                # Update particle velocity by averaging
                # current and new forces
                pp.leap_velocity(dt, 0.5*(force+force_new))
                p.leap_velocity(dt, 0.5*(force+force_new))

                # Re-define force value
                force = force_new


                #define new two particles separation
                pos = np.linalg.norm(pp.position)
                #define new two particles separation
                sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))

                KE += Particle.kin_energy(pp)
                PE += Particle.pot_energy(pp,p)/2
                energy = KE + PE

                # Increase time
                time += dt

                def VMD(no_particle, time, pp):                                         #ammend this later???
                    return str(no_particle) + '\n' + str(time) + '\n' + str(pp)         #cant call VMD function here??
                VMD = VMD(no_particle, time, pp)
                #print(VMD)

                # Output particle information
                trajectory_file.write(VMD)

                # Append information to data lists
                time_list.append(time)
                pos_list.append(pos)
                sep_list.append(sep_pos)
                energy_list.append(energy)


            #print(time_list)
            #Calculate Observables:                                                           #printing data list is weird here??
            # maxima_oxyeng = ext(np.array(energy_list), np.greater)[0]
            # print(maxima_oxyeng)






    # Post-simulation:
    # Close output file
    trajectory_file.close()

main()

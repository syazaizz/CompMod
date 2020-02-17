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

# def VMD(no_particle, time, pp):
#     return str(no_particle) + '\n' + str(time) + '\n' + str(pp)                             #check pp?? VMD format
    #return str("{0:d}\n".format(no_particle)) + str("{0}\n".format(time)) + str(pp)        #whats this?     return str(no_particle) + '\n'

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
    #force = []
    for i in range(0,no_particle):
         particle.append(Particle.from_file(particle_file))

    #make force list [0,..] for no particle
    force=np.zeros(no_particle)

    # for f in force:
    #     print(f)
    # for pp in particle:
    #     print(pp)
    # for ps in enumerate(particle):
    #     print(ps)

    no_part = str(no_particle)+'\n'
    step_no = str("step=0")+'\n'
    trajectory_file.write(no_part)
    trajectory_file.write(step_no)

    #COM correction
    M = 0
    P = 0
    for pp in particle:
        M += pp.mass
        P += pp.velocity*pp.mass
    v_com = P/M
    #print(v_com)

    KE = 0
    PE = 0

    pos_list = []
    sep_list = []                                             #might need to specify sep between sun, and moon earth?


    for i, pp in enumerate(particle):
        #print(i)
        #define new particle velocity
        pp.velocity -= v_com
        #print(pp.velocity)


        #write to trajectory_file
        #trajectory_file.write(vmd)
        part = str(pp)+'\n'
        trajectory_file.write(part)

        pos = np.linalg.norm(pp.position)
        pos_list.append(pos)
        #print(pos_list)

        sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))                            #how to exclude the 1st row?
        sep_list.append(sep_pos)


        # Write out initial energy and force
        KE += Particle.kin_energy(pp)                                                          #check units of G compatible


        for j, p in enumerate(particle):
            # print(j)
            if j > i:

                PE += Particle.pot_energy(pp,p)/2

                #for f in force:
                force += Particle.force(pp,p)    #cxheck unit G

    # print(KE)
    #print(KE)
    # for f in force:
    #     print(f)


    #initial total energy
    energy = KE + PE
    print(energy)



    #initialise data lists
    time_list = [time]
    energy_list = [energy]
    # print(pos_list)
    # print(sep_list)

    force_new = np.zeros(no_particle)

    for i in range(int(numstep)):
        no_part = str(no_particle).format(i)+'\n'
        step_no = "step={0}".format(i+1)+'\n'
        trajectory_file.write(no_part)
        trajectory_file.write(step_no)
        #print(energy)


        for j, pp in enumerate(particle):

            pp.leap_pos2nd(dt, force)
            pp.leap_velocity(dt, 0.5*(force+force_new))
            KE += Particle.kin_energy(pp)                                       #KE here too big???
            #print(KE)
            # part = str(pp)+'\n'
            # trajectory_file.write(part)

            #define new two particles separation
            pos = np.linalg.norm(pp.position)
            #define new two particles separation
            sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))


            for k, p in enumerate(particle):

                # Update particle position
                # pp.leap_pos2nd(dt, force)
                p.leap_pos2nd(dt, force)

                part = str(p)+'\n'
                trajectory_file.write(part)

                if k != j:

                    # Update force
                    force_new += Particle.force(pp,p)   #fx
                    PE += Particle.pot_energy(pp,p)/2
                #print(force_new)

                # Update particle velocity by averaging
                # current and new forces
                # pp.leap_velocity(dt, 0.5*(force+force_new))
                p.leap_velocity(dt, 0.5*(force+force_new))


                force = force_new

                # KE += Particle.kin_energy(pp)
                #if k > j:
                #if k != j:
                    #PE += Particle.pot_energy(pp,p)/2

        #print(KE)
        energy = KE + PE
        #print(energy)


        # Increase time

        time += dt


        # Append information to data lists
        time_list.append(time)
        #print(time)
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

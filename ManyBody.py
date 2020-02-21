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

def obs(list):
    max_pos = ext(np.array(list), np.greater)[0]       #apoapsis index
    min_pos = ext(np.array(list), np.less)[0]          #periapsis index
    max_val = list[max_pos[0]]                         #apoapsis value,       ## put [0] if more than one extrema ##
    min_val = list[min_pos[0]]                         #periapsis value
    t1 = time_list[max_pos[0]]                         # simulation need to have at least 2.5 complete orbit
    t2 = time_list[max_pos[1]]
    T = t2 - t1                                        # unit = day
    return (max_val, min_val, T)




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
    dt = np.array([float(i) for i in data[0:1]])   #data[0]
    numstep = np.array([float(i) for i in data[1:2]])
    time = 0.0
    print("dt =", dt, "and numstep =", numstep)


    #set up particles initial condition #list n loops to read p1, p2 ...
    particle_file = open(particle_file_name, "r")
    no_particle = count(particle_file)
    particle_file = open(particle_file_name, "r")
    particle=[]
    #no_particle = 4
    #force = []
    for i in range(0,no_particle):
         particle.append(Particle.from_file(particle_file))

    #make force list [0,..] for no particle

    force = np.zeros(3).reshape((3,))                                            ## check this force ???

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


    #pos_list = np.zeros(no_particle)
    #pos_list = []
    sep_list = []                                             ## might need to specify sep between sun, and moon earth ???
    sep_moon = []


    #for i, pp in enumerate(particle):
    #for pp in particle:
        #print(force)

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


        # pos = np.linalg.norm(pp.position)
        # pos_list.append(pos)
        #pos_list += pos

        pos = np.linalg.norm(pp.position)
        pos_list.append(pos)

        #print(pos_list)

        sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))                            #how to exclude the 1st row?
        sep_list.append(sep_pos)


        # eth_moon = np.linalg.norm(Particle.Vect_Sep(particle[3],particle[11]))                 # sep between earth n moon
        # sep_moon.append(eth_moon)


        # Write out initial energy and force
        KE += Particle.kin_energy(pp)                                                          #check units of G compatible


        for j, p in enumerate(particle):

        #for p in particle:
            p.velocity -= v_com
            #print(force)
            # if j != i:
            if j != i:
                # print(force)
                PE += Particle.pot_energy(pp,p)/2
                #print(PE)
                #for f in force:
                force += Particle.force(pp,p)    #cxheck unit G
                #print(force)

    # print(KE)
    # print(PE)
    # for f in force:
    #     print(f)
    # ps_list = np.array(pos_list)
    #print(sep_list)


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


    force_new = np.zeros(3).reshape((3,))

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

            part = str(pp)+'\n'
            trajectory_file.write(part)

            #define new particle position from center
            # pos = np.linalg.norm(pp.position)
            # pos_list.append(pos)


            #define new two particles separation
            sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))
            sep_list.append(sep_pos)

            eth_moon = np.linalg.norm(Particle.Vect_Sep(particle[3],particle[11]))                 # sep between earth n moon
            sep_moon.append(eth_moon)

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


                # part = str(p)+'\n'
                # trajectory_file.write(part)

                # if k > j:
                #
                #     new_pos = p.position - pp.position
                #     # Update force
                #     force_new += Particle.force(pp,p)   #fx
                #     #PE += Particle.pot_energy(pp,p)/2
                #     PE += Particle.pot_energy(pp,p)/2


                # part = str(p)+'\n'
                # trajectory_file.write(part)

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


        time += float(dt)
        time_list.append(time)
        #print(time)

        #print(KE)
        energy = KE + PE
        #print(energy)


        # Increase time


        #time += dt


        # Append information to data lists
        #time_list.append(time)
        #print(time)
        # pos_list.append(pos)
        # sep_list.append(sep_pos)

        time += dt


        # Append information to data lists
        time_list.append(time)
        #print(time)
        pos_list.append(pos)
        sep_list.append(sep_pos)
        energy_list.append(energy)






    #print(time_list)
    #print(pos_list)

    """
    #Calculate Observables
    # 1: Apo and Periapses
    # 2: Orbital Period

    sep = np.array([sep_list[i:i + no_particle] for i in range(0, len(sep_list), no_particle)])
    #print(sep)

    Sun = []
    Hg = []
    Ven = []
    Eth = []
    Mar = []
    Jup = []
    Sat = []
    Urn = []
    Nep = []
    Plu = []
    Hcom = []


    for i in range(0,len(sep)):
        lst = sep[i]

        Sun_pos = np.array([lst[0]])
        Hg_pos = np.array([lst[1]])
        Ven_pos = np.array([lst[2]])
        Eth_pos = np.array([lst[3]])
        Mar_pos = np.array([lst[4]])
        Jup_pos = np.array([lst[5]])
        Sat_pos = np.array([lst[6]])
        Urn_pos = np.array([lst[7]])
        Nep_pos = np.array([lst[8]])
        Plu_pos = np.array([lst[9]])
        Hcom_pos = np.array([lst[10]])

        Sun.append(Sun_pos)
        Hg.append(Hg_pos)
        Ven.append(Ven_pos)
        Eth.append(Eth_pos)
        Mar.append(Mar_pos)
        Jup.append(Jup_pos)
        Sat.append(Sat_pos)
        Urn.append(Urn_pos)
        Nep.append(Nep_pos)
        Plu.append(Plu_pos)
        Hcom.append(Hcom_pos)


    #print(Hg)
    apo_Hg, per_Hg, T_Hg = obs(Hg)
    apo_Ven, per_Ven, T_Ven = obs(Ven)
    apo_Eth, per_Eth, T_Eth = obs(Eth)
    apo_Mar, per_Mar, T_Mar = obs(Mar)
    apo_Jup, per_Jup, T_Jup = obs(Jup)
    apo_Sat, per_Sat, T_Sat = obs(Sat)
    apo_Urn, per_Urn, T_Urn = obs(Urn)
    apo_Nep, per_Nep, T_Nep = obs(Nep)
    apo_Plu, per_Plu, T_Plu = obs(Plu)
    apo_Hcom, per_Hcom, T_Hcom = obs(Hcom)


    # apo_moon, per_moon, T_moon = obs(sep_moon)


    # ps_list = np.array([pos_list[i:i + no_particle] for i in range(0, len(pos_list), no_particle)])
    #ps_list = np.array([pos_list[i] for i in range(0, len(pos_list))])
    # print(ps_list)


    # 3: Energy fluctuation

    # Plot energy to screen, to check the position of max energy
    pyplot.title('total energy vs time')
    pyplot.xlabel('Time, *10.18 fs')
    pyplot.ylabel('Energy, eV')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    maxima_eng = ext(np.array(energy_list), np.greater)[0]
    max = energy_list[maxima_eng[0]]                                   ## put [1] if 2nd maxima greater than first, check by plotting graph ##
    Eo = energy_list[0]
    delta_E = max - Eo
    E_fluct = abs(delta_E/Eo)
    print("Energy fluctuation = ", E_fluct)
    print("\n")

    """

        #print(time_list)
        #Calculate Observables:                                                           #printing data list is weird here??
        # maxima_oxyeng = ext(np.array(energy_list), np.greater)[0]
        # print(maxima_oxyeng)



    # Post-simulation:
    # Close output file
    trajectory_file.close()

main()

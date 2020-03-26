"""
Project Solar System: ParticleManyBody main module.
THe program uses Velocity Verlet time integration mehod to simulate the main bodies of our solar systems,
interacting through Newtonian gravity.

Version: 01/2020

Author: Syaza Izzati
"""

import sys
import math
import numpy as np
from Particle import Particle
from scipy.signal import argrelextrema as ext
import matplotlib.pyplot as pyplot
import linecache


# Count number of particle in particle.txt file
def count(particle_file):
    count = 0
    with particle_file as f:
        for line in f:
            count += 1
        return count


# Calculate observables
def obs(pos_list, time_list):

    max_pos = ext(np.array(pos_list), np.greater)[0]       # apoapsis index
    min_pos = ext(np.array(pos_list), np.less)[0]          # periapsis index
    max_val = pos_list[max_pos[0]]                         # apoapsis value: put [0] if more than one extrema, unit = km
    min_val = pos_list[min_pos[0]]                         # periapsis value, unit = km


    t1 = time_list[max_pos[0]]                             # simulation need to have at least 2.5 complete orbit
    t2 = time_list[max_pos[1]]

    T = t2 - t1                                            # unit = day

    print("apoapsis (km) = ", max_val)                     # print observables to terminal
    print("periapsis (km) = ", min_val)
    print("orbital period (day) = ", T)
    print("\n")
    return (max_val, min_val, T)


# Calculate observables for Neptune
def nept(pos_list, time_list):                             # function to calculate observables for neptune
    max_val = max(pos_list)
    min_val = min(pos_list)

    t1 = time_list[pos_list.index(max_val)]                # this method is more suitable to get accurate orbital period for neptune
    t2 = time_list[pos_list.index(min_val)]
    T = 2*abs(t2 - t1)

    print("apoapsis (km) = ", max_val)                     # print observables to terminal
    print("periapsis (km) = ", min_val)
    print("orbital period (day) = ", T)
    print("\n")
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
    energy_file = open("energy.dat", "w")


    # parameter of the simulation
    param_file = open(parameter_file, "r")                               # parameter_file contain data for parameter
    #line = param_file.readline()
    line = linecache.getline(parameter_file,1)                           # read specific no of line from parameter.txt
    data = line.split(",")
    dt = np.array([float(i) for i in data[0:1]])                         # dt in parameter_file, unit = day
    numstep = np.array([float(i) for i in data[1:2]])                    # number of steps for the simulation
    time = 0.0                                                           # simulation starts at time = 0 day
    print("dt (day) =", dt, "and numstep =", numstep)
    print("\n")


    # set up particles initial condition
    particle_file = open(particle_file_name, "r")
    no_particle = count(particle_file)                        # particle_file contains data for all particles: label, mass, position and velocity
    particle_file = open(particle_file_name, "r")

    particle=[]
    # no_particle = 12                                        # use this line to get certain number of particles from particle_file

    for i in range(0,no_particle):                            # create list of particles
         particle.append(Particle.from_file(particle_file))

    force = []
    for i in range(no_particle):                              # create list of force for each particle
        force.append(np.zeros(3))


    # write to trajectory file
    no_part = str(no_particle)+'\n'
    step_no = str("step=0")+'\n'
    trajectory_file.write(no_part)
    trajectory_file.write(step_no)


    # COM correction
    M = 0
    P = 0
    for pp in particle:
        M += pp.mass
        P += pp.velocity*pp.mass
    v_com = P/M


    KE = 0
    PE = 0
    sep_list = []
    sep_moon = []

    for i, pp in enumerate(particle):

        pp.velocity -= v_com                                                    # define new particle velocity

        part = str(pp)+'\n'                                                     # write to trajectory_file
        trajectory_file.write(part)

        sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))             # calculate separation between Sun and other planets inclusing halley comet
        sep_list.append(sep_pos)

        eth_moon = np.linalg.norm(Particle.Vect_Sep(particle[3],particle[11]))  # calculate separation between earth and moon
        sep_moon.append(eth_moon)

        KE += Particle.kin_energy(pp)                                           # calculate initial kinetic energy

        for j, p in enumerate(particle):

            if j != i:

                PE += Particle.pot_energy(pp,p)/2                               # calculate initial potential energy

                force_j_on_i = Particle.force(pp,p)                             # calculate force for each particle
                force[i] += force_j_on_i
                force[j] += -1*force_j_on_i                                     # Newton's 3rd law: avoid double counting

    for i in range(no_particle):
        force[i] = 0.5*force[i]


    # initial total energy
    energy = KE + PE

    # write to energy_file
    energy_file.write("Energy_data: (time [day], energy [kg*km^2/day^2])\n")
    energy_file.write("{0:f} {1:12.8f}\n".format(time,float(energy)))


    # initialise data lists
    time_list = [time]
    energy_list = [energy]


    #time integration
    for i in range(int(numstep)):

        if i%10==0:                                                             # only write to trajectory_file every 10 timesteps
            no_part = str(no_particle).format(i)+'\n'
            step_no = "step={0}".format(i+1)+'\n'
            trajectory_file.write(no_part)
            trajectory_file.write(step_no)

        KE = 0
        PE = 0

        # update position
        for j, pp in enumerate(particle):

            pp.leap_pos2nd(dt, force[j])

            if i%10==0:
                part = str(pp)+'\n'
                trajectory_file.write(part)

            # define new two particles separation
            sep_pos = np.linalg.norm(Particle.Vect_Sep(particle[0],pp))                     # separation between sun and other planets inclusing halley comet
            sep_list.append(sep_pos)

            eth_moon = np.linalg.norm(Particle.Vect_Sep(particle[3],particle[11]))          # separation between earth n moon
            sep_moon.append(eth_moon)

        # calculate force and potential energy on each particles
        force_new = []
        for j in range(no_particle):
            force_new.append(np.zeros(3))

        for j, pp in enumerate(particle):
            for k, p in enumerate(particle):
                if k != j:
                    PE += Particle.pot_energy(pp,p)/2

                    force_k_on_j = Particle.force(pp,p)
                    force_new[j] += force_k_on_j    #cxheck unit G
                    force_new[k] += -1*force_k_on_j #Newton's 3rd law: avoid double counting

        for i in range(no_particle):
            force_new[i] = 0.5*force_new[i]

        # update velocity and calculate kinetic energy
        for j, pp in enumerate(particle):
            pp.leap_velocity(dt, 0.5*(force[j]+force_new[j]))
            KE += Particle.kin_energy(pp)

        # re-define force value
        for j, pp in enumerate(particle):
            force[j] = force_new[j]

        # Increase time
        time += float(dt)
        time_list.append(time)

        # calculate total energy and write to energy_file
        energy = KE + PE
        energy_file.write("{0:f} {1:12.8f}\n".format(time,float(energy)))
        energy_list.append(energy)


    # Calculate Observables
    # 1: Apo and Periapses
    # 2: Orbital Period
    print("Observables")
    print("\n")


    # obtain data for separation between Sun and other planets including halley comet
    sep = np.array([sep_list[i:i + no_particle] for i in range(0, len(sep_list), no_particle)])
    list = [ [] for j in range(len(particle))]
    for i in range(1,len(sep)):
        lst = sep[i]
        dat = []
        for j in range(0, len(particle)):
            da = lst[j]
            dat.append(da)
            list[j].append(np.array([dat[j]]))
    for i in range(1, len(particle)-2):                                # exclude Moon and Neptune
        print (particle[i].label + ": ")
        obs(list[i], time_list)
    print("Neptune:")                                                  # calculate observables for Neptune
    nept(list[-2], time_list)


    # obtain data for Earth and Moon separation
    mn = np.array([sep_moon[i:i + no_particle] for i in range(0, len(sep_moon), no_particle)])
    mn_list = [ [] for j in range(len(particle))]
    for i in range(1,len(mn)):
        mn_lst = mn[i]
        mn_dat = []
        for j in range(0, len(particle)):
            mn_da = mn_lst[j]
            mn_dat.append(mn_da)
            mn_list[j].append(np.array([mn_dat[j]]))
    print("Moon:")
    obs(mn_list[0], time_list)                                         # will be the same for [1],[2],..,[len(particle)]



    # 3: Energy fluctuation

    # Plot energy to screen, to check the position of max energy and min energy
    pyplot.title('total energy vs time')
    pyplot.xlabel('Time, [day]')
    pyplot.ylabel('Energy, [kg*km^2/day^2]')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    minima_eng = ext(np.array(energy_list), np.less)[0]
    min = energy_list[minima_eng[0]]                                   # put [1] if 2nd minima greater than first: check by plotting graph
    Eo = energy_list[0]
    delta_E = abs(min - Eo)
    E_fluct = abs(delta_E/Eo)
    print("Initial total energy, Eo (kg*km^2/day^2) = ", energy)
    print("Energy fluctuation (kg*km^2/day^2) = ", E_fluct)
    print("\n")


    # Post-simulation:
    # Close output file
    trajectory_file.close()
    energy_file.close()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

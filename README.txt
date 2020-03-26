==============================
Project Solar System: Simulation of N-body systems interacting through Newtonian gravity for the main bodies of our solar system,
starting from realistic initialconditions.

Author: Syaza Izzati
Matriculation no: s1745516
Version: 03/2020

==============================
RUNNING THE CODE:

In the correct directory, type to terminal:
python3 ManyBody.py parameter.txt planet.txt traj.xyz

Input file:
a) parameter.txt: contains numerical parameters of the simulation: timesteps, dt (day) and number of steps, numstep.
  - dt = 1 day, numstep = 181140. The value of numstep was decided based on the orbital period of Pluto (outer planet in solar system),
    the simulation need to have at least 2.5 complete orbit for observables calculation.
b) planet.txt: contains the details of the particles: labels, masses(kg), starting positions(km), starting velocities(km/day).

Output file:
a) traj.xyz: contains string representing the particle in XYZ-compatible format for a VMD trajectory.
b) energy.dat: contains the list of total energy with time, to track the total energy of the system.

==============================
1. Particle.py:
Contains class that describe the particles moving in 3D, reading the particles data from "planet.txt" file.

Unit:
mass = kg
position = km
velocity = km/day
timestep, dt = day
Gravitational constant, G = 6.674 30 x 10-11 m^3 kg^(-1) s^(-2)

==============================
2. ManyBody.py:
A module that contains Velocity Verlet time integration that describe N-body systems interacting through Newtonian gravity,
and is used to simulate the main bodies of our solar system.

Output:
a) Prints the value of timestep, dt and number of steps, numstep to terminal.
b) Produces and prints the observable of each planet including Halley Comet and Moon to terminal: the apoapsis, periapsis and orbital period.
c) Produces plots of total energy energy as function of time. Also save the energy list to a file: energy.dat
d) Prints the total energy, Eo and energy fluctuation of the system to terminal.

==============================

"""
 Project Solar System: Particle, a class to describe 3D planets

 Version: 01/2020

 Author: Syaza Izzati
"""
import numpy as np

class Particle(object):
    """
    Class to describe 3D planets.

    Properties:
    position(numpy arrays) - position along the x, y and z axis
    velocity(numpy arrays) - velocity along the x, y and z axis
    mass(float) - particle mass
    label(string) - particle label

    Methods:
    * formatted output
    * first-order velocity update
    * second order position updates
    * kinetic energy
    """

    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle3D instance

        :param pos: position as numpy array (unit - km)
        :param vel: velocity as numpy array (unit - km/day)
        :param mass: mass as float(file_handle) (unit - kg)
        :param label: label as string
        """
        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label



    def __str__(self):
        """
        Define output format.
        For particle p=([x_pos, y_pos, z_pos], [x_vel, y_vel, z_vel], mass, label) this will print as
        "label x_pos y_pos z_pos"
        """
        return str(self.label) +" "+ str(self.position[0]) +" "+ str(self.position[1]) + " " + str(self.position[2])


    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*f(t)/m

        :param dt: timestep as float
        :param force: force on particle as float
        """
        self.velocity +=  dt*force/self.mass


    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 0.5*dt^2*f(t)/m

        :param dt: timestep as float
        :param force: current force as float
        """
        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass


    #kinetic energy
    def kin_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return 0.5*self.mass* (np.linalg.norm(self.velocity))**2


    #read content from filename
    # call Particle __init__ method
    @staticmethod
    def from_file(file_handle):

        line = file_handle.readline()                            # read line in file_handle
        data = line.split(",")                                   # data in file_handle is separated by symbol: ","
        label = data[0].strip('\n')                              # label of particle in file_handle

        data = [float(data_val) for data_val in data[1:9]]
        if len(data) == 0:
            return 0

        else:
            pos = np.array([float(i) for i in data[1:4]])        # position of particle in file_handle
            vel = np.array([float(i) for i in data[4:7]])        # velocity of particle in file_handle
            mass = np.array([float(i) for i in data[0:1]])       # mass of particle in file_handle
            return Particle(pos, vel, mass, label)


    #relative vector separation of two particles
    @staticmethod
    def Vect_Sep(p1,p2):
        return p1.position - p2.position

    #gravitational potential energy between two particles
    @staticmethod
    def pot_energy(p1,p2):
        """
        Return potential energy as
        -G*mass_1*mass_2/|r_12|

        G: gravitational constant, 6.674 30 x 10-11 m^3 kg^(-1) s^(-2)
        r_12: vector separation between particle 1 and particle 2
        """
        return -(4.982338253E-10)*p1.mass*p2.mass/np.linalg.norm(Particle.Vect_Sep(p1,p2))            # convert unit of G to km^3 kg^(-1) day^(-2) to match with the input units from file_handle

    #gravitational force between two particles
    @staticmethod
    def force(p1,p2):
        """
        Return force as
        (-G*mass_1*mass_2/|r_12|^3)*(r_12)

        G: gravitational constant, 6.674 30 x 10-11 m^3 kg^(-1) s^(-2)
        r_12: vector separation between particle 1 and particle 2
        """
        return (-(4.982338253E-10)*p1.mass*p2.mass/(np.linalg.norm(Particle.Vect_Sep(p1,p2)))**3)*Particle.Vect_Sep(p1,p2)         # convert unit of G to km^3 kg^(-1) day^(-2) to match with the input units from file_handle

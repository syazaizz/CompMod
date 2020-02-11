"""
 Project Solar System: Particle, a class to describe 3D planets

 Version: 01/2020
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

        :param pos: position as numpy array
        :param vel: velocity as numpy array
        :param mass: mass as float(file_handle)
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
        #return "label = " + str(self.label) +" "+ "x_pos = " + str(self.position[0]) +" "+"y_pos = " + str(self.position[1]) + " " +"z_pos = " + str(self.position[2])

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

        line = file_handle.readline()
        data = line.split(",")
        label = data[-1].strip('\n')

        #print(data)

        data = [float(data_val) for data_val in data[0:-1]]
        if len(data) == 0:
            return 0  #false
        # print(data)
        else:
            pos = np.array([float(i) for i in data[0:3]])
            vel = np.array([float(i) for i in data[3:6]])
            #print(vel)
            mass = np.array([float(i) for i in data[6:7]])
            #print(mass)
            #print(label)
            return Particle(pos, vel, mass, label)


    #relative vector separation of two particles
    @staticmethod
    def Vect_Sep(p1,p2):
        return p1.position - p2.position

    #gravitational potential energy between two particles
    @staticmethod
    def pot_energy(p1,p2):
        return -1*p1.mass*p2.mass/np.linalg.norm(Particle.Vect_Sep(p1,p2))

    #gravitational force between two particles
    @staticmethod
    def force(p1,p2):
        return (-1*p1.mass*p2.mass/(np.linalg.norm(Particle.Vect_Sep(p1,p2)))**3)*Particle.Vect_Sep(p1,p2)

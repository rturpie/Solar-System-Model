#  ------ Computer Modelling Excersise 2 ---------

# importing needed libraries

import numpy as np

# creating a Python class Particle3D that describes a point-like particle moving in 3D space

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    label(string) - particle's label
    position(array) - position along the x,y and z axis
    velocity(array) - velocity along the x, y and z axis
    mass(float) - particle mass

    Methods:
    * initialise
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    * create a particle from a file handle argument
    * calculate particle seperation
    """
    def __init__(self, label, pos, vel, mass):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param pos: position as an array
        :param vel: velocity as an array
        :param mass: mass as float
        """

        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label

    def __str__(self):
        """
        Define output format.
        For particle with position (2.0, 0.5, 1.0) this will return the string
        "2.0 0.5 1.0"
        """

        # Note the divisions by 10**9 are to stop VMD from struggling to deal with huge numbers that it is not used to. It will still work. 

        return str(self.label) + " " + str((self.position)[0]/10**9) + "  "  + str((self.position)[1]/10**9) + " " + str((self.position)[2]/10**9)

    def kinetic_energy(self):
        """
        Return kinetic energy as float
        """
        velocity_squared = np.dot(self.velocity,self.velocity)

        return 0.5*self.mass*velocity_squared

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle vector as numpy array
        """

        self.velocity = self.velocity + dt*force/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """

        self.position = self.position + dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force vector as numpy array
        """

        self.position = self.position + dt*self.velocity + 0.5*dt**2*force/self.mass

    @staticmethod
    def create_particle(file_handle):
        """
        create a particle from a file entry

        :param file_handle: the file handle the paramenters of the particle will be taken from
        :return: the new particle as a Particle 3D object
        """

        filein = open(file_handle,"r")
        lines = filein.readlines()
        filein.close()

        system = []

        for l in lines:
            l = l.strip()
            details = l.split(",")
            l = details[0]
            p = np.array([float(details[1]),float(details[2]),float(details[3])])
            v = np.array([float(details[4]),float(details[5]),float(details[6])])
            m = float(details[7])
            particle = Particle3D(l,p,v,m)

            system.append(particle)

        return system


    @staticmethod
    def particle_separation(particle1, particle2):
        """
        calculate the relative vector seperation of two particles

        :param particle1: the first particle
        :param particle2: the second particle
        :return: the relative vector seperation of the first and second particle
        """
        return particle1.position - particle2.position

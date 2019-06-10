import numpy as np
from Particle3D import Particle3D


def Find_potential_energy(particle1, particle2):

    """
    Method to return the gravitational potential energy between two objects

    Potential energy is given by Ep = -Gm1m2/r

    :param particle1: Particle3D instance
    :param partcile2: Particle 3D insatnce
    :return: The potential energy between two particles
    """

    r12 = Particle3D.particle_separation(particle1, particle2)
    rscalar = np.linalg.norm(r12)

    G = 6.67408E-11        # Gravitational constant

    potential_energy = (-G*particle1.mass*particle2.mass)/rscalar

    return potential_energy



def Find_kinetic_energy(particle1):

    """
    Method to return the kinetic energy of an object,
    Kinetic energy is given by Ek = 0.5mv^2

    :param particle1: Particle3D instance
    :return: The kinetic energy of the particle
    """

    kinetic_energy = 0.5*particle1.mass*(np.linalg.norm(particle1.velocity)**2)

    return kinetic_energy


def Find_force_vector(particle1, particle2):

    """
    Method to return the force due to gravity as a vector between two objects.

    Force is given by F = (Gm1m2/r^2)rhat

    Where rhat is the unit vector between the two objects.

    :param particle1: particle3D instance
    :param partcile2: particle 3D instance
    :return: force between particle1 and particle2 as a numpy array
    """

    r12 = Particle3D.particle_separation(particle1, particle2)  # displacement between the particles, a vector
    rscalar = np.linalg.norm(r12)                               # distance between the particles, a scalar

    rhat = r12/rscalar      # unit vector
    
    G = 6.67408E-11         # Gravitational constant           

    force = -((G*particle1.mass*particle2.mass)/(rscalar**2))*rhat

    return force


def Find_all_forces(lst):

    """
    A method to calculate the forces on each object due to every other object using the Find_Force_Vector function.

    :param lst: a list of Particle3D objects
    :return: a list of numpy arrays, where each array is the total force acting on the body of that indexing 
    """

    all_total_forces = []       # forces acting on each particle
    all_individual_forces = []  # will be a list of lists of the forces acting on a particle from each other particle, or a zero if not needed


    for i in range(0, len(lst)):
        individualforces = []   # list to go into list of lists
        totalforce = np.array([0,0,0])

        for j in range(0,len(lst)):
            # ie particle and itself
            if j == i : 
                individualforces.append(0)

            # Calculate force between the particles
            elif j > i:
                added = Find_force_vector(lst[i],lst[j])
                totalforce = totalforce + added
                individualforces.append(added)

            # use Newton's third law as force has already been calculated for the other particle
            elif j < i:
                totalforce = totalforce - all_individual_forces[j][i]
                individualforces.append(0)

        all_individual_forces.append(individualforces)

        all_total_forces.append(totalforce)


    return all_total_forces


def Update_velocity(lst,force_lst,dt):

    """
    A  method to update the velocities of particles stored in a list.

    :param lst: a list containing Particle3D instances
    :param force_lst: a list of the forces acting on each object
    :param dt: timestep as float
    :return: the updated list of velocities which are numpy arrays
    """

    updated_lst = []

    for i in range(0, len(lst)):

        lst[i].leap_velocity(dt,force_lst[i])   #calculating updated velocity of particle

        updated_lst.append(lst[i])

    return updated_lst


def Update_position(lst,force_lst,dt):

    """
    A method to update the position of particles in a list.

    :param lst: a list containing Particle3D instances
    :param force_lst: a list of the forces acting on each object
    :param dt: timestep as float
    :return: the updated list of positions which are numpy arrays
    """
    updated_position_lst = []

    for i in range(0, len(lst)):

        updated_position_lst.append(lst[i].leap_pos2nd(dt,force_lst[i]))

    return updated_position_lst


def Calculate_total_energy(lst):

    """
    A method to calculate the total energy of the system due to particle pair interactions (potential energies) and their kinetic energies,
    using the functions Find_kinetic_energy and Find_potential_energy.

    :param lst: a list containing Particle3D instances representing the system
    :return: 3-tuple of floats, (total energy, kinetic energy, potential energy)
    """

    total_energy = 0
    potential_energy = 0
    kinetic_energy = 0

    for l in lst:

        add_kin = Find_kinetic_energy(l)
        kinetic_energy += add_kin
        total_energy += add_kin # adds kinetic energy of each particle3D to total energy

        for i in range(lst.index(l) + 1, len(lst)):

            add_pot = Find_potential_energy(l,lst[i])
            potential_energy += add_pot
            total_energy += add_pot # adds potential between that particle and each further particle in lst to total energy

    return total_energy,kinetic_energy,potential_energy

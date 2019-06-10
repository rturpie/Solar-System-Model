# import needed libraries and functions from other modules

from Particle3D import Particle3D
import numpy as np
import sys
import calculation as calc
import VMD
import matplotlib.pyplot as pyplot
import math


def Apsides(distance_list):

    """
    Method to return the apo- and periapses of two bodies
    apoapsis is given by the longest distance between the orbitting bodies during rotation
    periapsis is given by the shortest distance between the orbitting bodies during rotation

    :param distance_list: list of scalar distances between the two bodies
    :return: the apoapsis and periapsis of the orbit as floats, respectivly
    """

    return max(distance_list), min(distance_list)


def orbit(satellite_list, parent_body_list, time_list):

    """
    Method to return the orbital period of each body around its parent body.
    The total angle it sweeps out is calculated by adding the angle swept out between successive positions at successive times in time_list.
    The total time that has passed is then divided by the number of rotations (ie total angle swept / 2pi) to give the period
    If one full orbit has not been completed, a string is returned which later is used to prompt the user to give more time steps,
    to avoid innacurate orbital periods

    :param satellite_list: list of positions of orbitting body as numpy arrays
    :param parent_body_list: list of positions of orbitted body as numpy arrays
    :param time_list: list of times that correspond to each position of the two bodies
    :return: the orbital period as a float, otherwise a string if one orbit has not been completed in the simulation
    """   

    # initialise and the create list of the vector difference between the two bodies as numpy arrays

    vector_diff_list = []

    for i in range(len(satellite_list)):
        vector_diff_list.append(satellite_list[i] - parent_body_list[i])

    orbitals_time = 0
    total_angle = 0

    for j in range(len(vector_diff_list) - 1):
        a = vector_diff_list[j]  # vector difference in position at this current point
        b = vector_diff_list[j+1] # vector difference in position at next point
        theta = math.acos((np.dot(a,b)/((np.linalg.norm(a))*(np.linalg.norm(b))))) # angle between a and b
        total_angle += theta
        orbitals_time = 10*time_list[j]

    if total_angle < 2*math.pi:
        return "try more time steps"

    else:

        return (orbitals_time/(total_angle/2*math.pi))/(60*60*24)



def add_to_pos_lists(lst,d,p):

    """
    Method to update lists in d and p at different points of time in the simulation, 
    where d stores within it the list of distances between each body in lst(except the sun) and its parent body,
    and p stores within it the list of positions as numpy arrays of each body in lst.
    Note that these lists being updated are the second element of pairs in the lists of pairs, d and p, where the first is the label of the body that the information is relevant to.
    Ie d and p are lists of 2-tuples of a string and a list of floating point numbers.

    :param lst: list of the bodies as Particle3D instances
    :param d: list of pairs, in which the first pair element is the body label as a string, and the second is a list of all the scalar distances of the body from the body it is orbitting
    :param p: list of pairs, in which the first pair element is the body label as a string and the second is a list of positions of the body as a numpy array
    """   

    # find the two bodies that are being orbitted and assign them to Sun and Earth
    # done in seperate loop so that the order of particles in the input file doesn't matter

    for l in lst:
        if l.label == 'Sun': 
            Sun = l
        if l.label == 'Earth':
            Earth = l

    for l in lst:

        p[lst.index(l)][1].append(l.position) # appends position of particle l to corresponding list of its positions in p
        
        if l.label == 'Sun': # sun is not orbitting anything so can ignore distances from anything else in this case
            pass
        elif l.label == 'Moon': # ensure distance from moon to earth is appended to the correct list in d  
            d[lst.index(l)][1].append(np.linalg.norm(Particle3D.particle_separation(l,Earth)))

        else: # everything else is orbitting sun so the distance between sun and that body is appended do the correct list in d
            d[lst.index(l)][1].append(np.linalg.norm(Particle3D.particle_separation(l,Sun)))


def com_correction(lst):

    """
    Method to correct the centre-of-mass motion for the n-body system
    Each body in the system has the centre of mass velocity (Σmivi)/(Σmi) taken from its initial velocity
    mi are all the masses of the system and vi all the initial velocities

    :param lst: a list containing Particle3D instances
    :return: a list containing the Particle3D instances with their corrected velocities
    """

    total_momentum = 0
    total_mass = 0
    corrected_system = []

    for l in lst:
        total_momentum += l.mass*l.velocity
        total_mass += l.mass

    v_com = total_momentum/total_mass

    for l in lst:
        l.velocity = l.velocity - v_com
        corrected_system.append(l)

    return corrected_system

# main method

def main():

    # Read name of input and output files from command line
    if len(sys.argv)!= 6:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0]  + " <input file of simulation parameters>" + " <input file of particle details>" + " <output vmd file>" + " <orbital outputfile>" + "<energy output file>")
        quit()
    else:
        
        parameter_file = sys.argv[1]
        particle_file = sys.argv[2]
        outfile_name = sys.argv[3]
        orb_file = sys.argv[4]
        energy_file = sys.argv[5]

    # create bodies
    # the input file has all positions in metres and velocities in ms^-1

    system = Particle3D.create_particle(particle_file)
    
    # centre of mass correction for the bodies

    com_correction(system)

    # begginings of position and distance lists for calculating orbitals and apsis

    positions = []
    distances = []

    for s in system:
        positions.append((s.label,[]))
        distances.append((s.label,[]))

    add_to_pos_lists(system,distances,positions)

    # read in simulation parameters from parameter file

    parameters = open(parameter_file, "r")
    lines = parameters.readlines()
    parameters.close()

    params = []

    for l in lines:
        l = l.strip()
        param, value = l.split(":")
        params.append((param,value))

    dt = list(map(lambda x: int(x[1]),(filter(lambda x: "dt" == x[0], params))))[0] # timestep in seconds
    n = list(map(lambda x: int(x[1]),(filter(lambda x: "n" == x[0], params))))[0] # for every n timesteps we will write to vmd, calculate energies etc (no unit)
    num_steps = list(map(lambda x: int(x[1]),(filter(lambda x: "num_steps" == x[0], params))))[0] # number of timesteps (no unit)

    # open vmd file to write to and initialise that it is on its first point

    vmd_file = open(outfile_name , "w")
    point = 1

    # initialise time

    time = 0

    # get initial forces 

    forces = calc.Find_all_forces(system)

    # get inital energies
    
    energies = calc.Calculate_total_energy(system)
    tot_energy = energies[0]
    kin_energy = energies[1]
    pot_energy = energies[2]

    # append energies and time to corresponding lists

    time_list = [time]
    tot_energy_list = [tot_energy]
    kin_energy_list = [kin_energy]
    pot_energy_list = [pot_energy]

    # start the time integration loop

    for i in range(num_steps):

        # update particle position
        calc.Update_position(system,forces, dt)

        # update forces
        new_forces = calc.Find_all_forces(system)

        # update particle velocities by averaging current and new forces
        
        average_forces = []
        
        for i in range(0,len(forces)):
            average_force = (forces[i] + new_forces[i])/2
            average_forces.append(average_force)

        calc.Update_velocity(system,average_forces, dt)

        # re define force
        forces = new_forces

        # increse time
        time += dt

        # every n timesteps

        if time % n == 0:

            VMD.Write_vmd(system,point,vmd_file) # write to vmd file
            point += 1 # update point
            
            # get energies

            energies = calc.Calculate_total_energy(system)
            tot_energy = energies[0]
            kin_energy = energies[1]
            pot_energy = energies[2]
            
            #update energy and time lists

            tot_energy_list.append(tot_energy)
            kin_energy_list.append(kin_energy)
            pot_energy_list.append(pot_energy)
            time_list.append(time)
            
            # append to position and distance lists

            add_to_pos_lists(system,distances,positions)

    # close vmd file
    
    vmd_file.close()

    # open orbital output file and write to it the apsides and orbital periods

    orb_file = open(orb_file , "w")

    for p in positions:
        
        if p[0] == 'Moon': # ensure moon orbits Earth
            
            orbitting = 0
            
            for i in positions:
                if i[0] == 'Earth':
                    orbitting = i

            # check at least one orbit has occured and if not suggests larger number steps

            test = orbit(p[1],orbitting[1],time_list)

            if isinstance(test, str):
                orb_file.write("For " + str(p[0]) + " try a larger number of steps to allow completion of at least one orbit. \n")

            # otherwise prints relevant inormation to orbital file

            else: 
                orb_file.write("The orbit of the " + str(p[0]) + " is: ")
                orb_file.write(str(orbit(p[1],orbitting[1],time_list)) + " days, \n")
                orb_file.write("and the apo- and periapses are (respectively): ")
                orb_file.write(str(Apsides(distances[positions.index(p)][1]))+ " metres. \n")

        elif p[0] == 'Sun': # ignores sun 
            pass
        
        else: # makes all other calculations on assumption of orbitting sun
            for i in positions:
                if i[0] == 'Sun':
                    orbitting = i

            # check at least one orbit has occured and if not suggest larger number steps

            test = orbit(p[1],orbitting[1],time_list)

            if isinstance(test, str): 
                orb_file.write("For " + str(p[0]) + " try a larger number of steps to allow completion of at least one orbit. \n")

            # otherwise prints relevant inormation to orbital file

            else: 
                orb_file.write("The orbit of " + str(p[0]) + " is: ")
                orb_file.write(str(orbit(p[1],orbitting[1],time_list)) + " days, \n")
                orb_file.write("and the apo- and periapses are (respectively): ")
                orb_file.write(str(Apsides(distances[positions.index(p)][1]))+ " metres. \n")
            
    # close orbital output file

    orb_file.close()

    energy_file = open(energy_file , "w")

    energy_file.write("The total energies, kinetic energies, potential energies (in Joules) are as follows: \n")
    for i in range(len(tot_energy_list)):
        energy_file.write(str(tot_energy_list[i]) + "," + str(kin_energy_list[i]) + "," + str(pot_energy_list[i]) + "\n")

    # close energy file

    energy_file.close()
    

    # plotting total energy, kinetic energy and potential energy of the system against time

    pyplot.title('Energy Fluctuation Solar System Total Energy vs Time')
    pyplot.xlabel('Time,s')
    pyplot.ylabel('Total Energy, J')
    pyplot.plot(time_list, tot_energy_list)
    pyplot.show()

    pyplot.title('Energy Fluctuation Solar System Kinetic Energy vs Time')
    pyplot.xlabel('Time,s')
    pyplot.ylabel('Kinetic Energy, J')
    pyplot.plot(time_list, kin_energy_list)
    pyplot.show()

    pyplot.title('Energy Fluctuation Solar System Potential Energy vs Time')
    pyplot.xlabel('Time,s')
    pyplot.ylabel('Potential Energy, J')
    pyplot.plot(time_list, pot_energy_list)
    pyplot.show()

# run the main

main()

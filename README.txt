
Instructions for use of ParticleManyBody.py



The files provided should be ordered in the terminal as follows:

python3 ParticleManyBody.py param.txt create.txt <some_file_name.xyz> <outfile_1.txt> <outfile_2.txt>

param.txt is the input file containing the simulation parameters. The top number, dt, is the timestep. The second number n dictates how often the updated particle information is written to the VMD file. The third number, numstep, is the number of steps that will be calculated in the time integration loop.

create.txt is the input file containing the initial conditions for all the bodies. The information for each body is contained within a single line. They are ordered from left to right: label, x position, y position, z position, x velocity, y velocity and finally z velocity.

some_file_name.xyz is the output file that the information for VMD will be written to.

outfile_1.txt is the output file that will contain the orbital periods and the apsides of the bodies calculated by the programme.

outfile_2.txt is the file that will contain the total energies, kinetic energies and potential energies of the system at various timesteps.

Included also is animate.xyz,data.txt and energies.txt, files we calculated using param.txt and create.txt, just as an example. Also feel free to change the parameters in param.txt and have fun:)

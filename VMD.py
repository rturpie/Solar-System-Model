
from Particle3D import Particle3D
import numpy as np
import calculation


def Write_vmd(lst,pnt,outfile):

    """
    Method to write to an output file information in the desired VMD format which will then be plotted using VMD to visualise the simulation.
    I.e. in the following format for each entry:
    
    n
    Point = 1
    s1 x11  y11  z11
    :  :    :    :
    sn xn1  yn1  zn1


    :param lst: bodies as list of Particle3D instances
    :param pnt: the point of this VMD entry
    :param outfile: the filehandle of the file the VMD information will be written to
    :return: prints entry for a timestep to the VMD file
    """

    outfile.write(str(len(lst)))
    outfile.write("\n")
    outfile.write("Point = " + str(pnt))
    outfile.write("\n")
    
    for l in lst:

        outfile.write(l.__str__()) # using Particle3D string method
        outfile.write("\n")


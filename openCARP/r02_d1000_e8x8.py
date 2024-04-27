#!/usr/bin/env python

import os # This module provides a portable way of using operating system dependent functionality
EXAMPLE_DESCRIPTIVE_NAME = 'Example01'
EXAMPLE_AUTHOR = 'Pascual'
EXAMPLE_DIR = os.path.dirname(__file__) # Return the directory name of pathname path
GUIinclude = True

from datetime import date # The datetime module supplies classes for manipulating dates and times

# carputils is a Python package containing several openCARP-related functionalities and a framework for running electrophysiology simulations and testing their outputs en masse
from carputils import settings
from carputils import tools
from carputils import mesh
from carputils.carpio import txt

import numpy as np # numpy is library for the Python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays

def parser():
    parser = tools.standard_parser() # Generate a standard argument parser for collection of common options
    group  = parser.add_argument_group('Experiment specific options') # The add_argument_group() method returns an argument group object which has an add_argument() method just like a regular ArgumentParser. When an argument is added to the group, the parser treats it just like a normal argument, but displays the argument in a separate group for help messages. The add_argument_group() method accepts title and description arguments
    
    # Filling an ArgumentParser with information about program arguments is done by making calls to the add_argument() method
    group.add_argument('--duration',
                        type = float,
                        default = 1000.,
                        help = 'Duration of simulation in [ms] (default: 1000.)')
    group.add_argument('--S1-strength',
                        type = float,
                        default = 200.,
                        help = 'pick transmembrane current stimulus strength in [uA/cm^2] = [pA/pF] considering fixed Cm (default: 20.)')
    group.add_argument('--S1-dur',
                        type = float,
                        default = 2.,
                        help = 'pick transmembrane current stimulus duration in [ms] (default: 2.)')
    return parser

def jobID(args):
    # Generate name of top level output directory
    today = date.today()
    return '{}_basic_{}'.format(today.isoformat(), args.duration)

@tools.carpexample(parser, jobID)
def run(args, job):
   
    # Block with a size of 50x50x1 mm and a spatial resolution of 0.2 mm
    # It is created: ((50/0.2)+1)*((50/0.2)+1)*((1/0.2)+1) = 378006 points that modelate the tissue mesh
    geom = mesh.Block(size=(50, 50, 1), resolution=0.2) # geom is the name of the variable that contains the absolute path where it is stored a file of the mesh. The inputs are the tissue size and the resolution
    geom.set_fibres(0, 0, 90, 90) # Set endocardial and epicardial fibre angle to 0 and endocardial and epicardial sheet angle to 90
    
    # Describe an axis-aligned cuboid for mesh tag assignment
    reg = mesh.BoxRegion((-25,-25,-0.5), (25,25,0.5), tag=1)
    geom.add_region(reg)

    # Generate the tissue mesh
    meshname = mesh.generate(geom)

    # Creation of the electrodes at a certan distance from the tissue mesh. It must be taken into account that each electrode spans a region of the tissue mesh created before.
    # It is created: ((21/3)+1)*((21/3)+1) = 64 electrodes that simulates the 8x8 electrode distribution of this simulation
    points = mesh.Block(size=(21, 21, 0), resolution=3) #It is created an absolute path where it is stored a file of the mesh 
    points.set_fibres(0, 0, 0, 0) # Set endocardial and epicardial fibre angle to 0 and endocardial and epicardial sheet angle to 0
    
    # Describe an axis-aligned cuboid for mesh tag assignment
    regelec = mesh.BoxRegion((-11,-11,-0.1), (11,11,0.1), tag=0)
    points.add_region(regelec)
    
    # Generate the electrodes distribution
    meshname2 = mesh.generate(points)
    
    # The geometry of the mesh to obtain the signal sources is centered in (0,0,0) and also the electrodes distribution but it is wanted to move the position in the z axis on the surface of the mesh called "geom" that modelates the tissue because this simulation is non-contact catheter mapping. In this simulation, it is wanted to have the electrode distribution at a distance of 0.5 mm from the tissue mesh
    # As the tissue mesh has a width of 1 mm, and it is located 0.5 mm on the positive axis of z and 0.5 on the negative axis of z, then 1 mm must be added to the third dimension of the electrodes distribution so that they are located at a 0.5 mm distance over the tissue mesh
    dir_name = points.mesher_opts(meshname2)[25] # It is obtained the name of the directory where it is stored the file that has the location of each electrode
    points_text = np.loadtxt(dir_name + '.pts',skiprows = 1) #It is loaded the file where it is stored the geometry of the source points (.pts file)
    
    points_text[:,2] = points_text[:,2] + 1000 # It is added the position of the 3rd column of the file to 1000 µm = 1 mm

    # Save the file of the mesh points fixed
    print(dir_name +'.pts')
    f=open(dir_name +'.pts','w')
    f.write('{}\n'.format(points_text.shape[0])) #openCarp
    np.savetxt(f, points_text, delimiter=' ')
    f.close()
    
    # Electrograms setup
    # phie_rec_ptf = dir_name --> Defines the basename for the phie recovery file
    # phie_recovery_file = electrograms --> Defines the file name used to output recovered extracellular potentials
    # phie_rec_meth = 1 --> Defines the method for recovering phie with monodomain runs
    egm = ['-phie_rec_ptf', dir_name ,
    '-phie_recovery_file', 'electrograms',
    '-phie_rec_meth', 1]

    # Stimulus options: one stimulus called S1 of transmembrane current (type 0) with a strength and duration established by the user (if the user does not set any value, the default values are 200 uA/cm^2 and 2 ms)
    stim = ['-num_stim',             1,
            '-stimulus[0].name',    'S1',
            '-stimulus[0].stimtype', 0, 
            '-stimulus[0].strength', args.S1_strength,
            '-stimulus[0].duration', args.S1_dur ]

    # Generate a boundary condition definition
    # def carputils.mesh.block.block_boundary_condition(block, entity, index, coord, lower = True, bath = False, verbose = False)
        #block	Block The geometry to generate boundary conditions for
        #entity	str The openCARP boundary condition type to use (e.g. 'stimulus')
        #index	int The boundary condition index
        #coord	str The direction perpendicular to which the BC is generated. One of ('x', 'y', 'z').
        #lower	bool, optional True if to generate BC at lower bound of coord (default), False for upper bound
        #bath	bool, optional Generate BC at bounds of tissue+bath if True, tissue only otherwise	
    electrode = mesh.block_boundary_condition(geom, 'stimulus', 0, 'x', True)

    # Ionic setup (Coutermanche model)
    imp = ['-num_imp_regions',1, '-imp_region[0].im','COURTEMANCHE','-imp_region[0].num_IDs',1, '-imp_region[0].ID',1]
    
    # Conduction velocity is assigned in different directions (longitudinal, transversal and normal) for the component of the tissue 
    cv=['-num_gregions',1,'-gregion[0].num_IDs',1,'-gregion[0].ID',1,'-gregion[0].g_il',0.1438,'-gregion[0].g_it',0.1438,'-gregion[0].g_in',0.1438]

    # Get basic command line, including solver options (ionic setup can be in other file)
    cmd = tools.carp_cmd() #Constructs a list of command line arguments. This will automatically include the loading of the correct options for the specified solvers

    # Simulation setup
    cmd += ['-simID',    job.ID,
            '-meshname', meshname,
            '-dt',       50,
            '-tend',     args.duration]

    # Adding al the features
    cmd += stim + electrode + imp + cv + egm

    if args.visualize:
        cmd += ['-gridout_i', 3]
        cmd += ['-gridout_e', 3]
        cmd += ['-spacedt', 0.1]


    # Run simulation 
    job.carp(cmd) 

    
if __name__ == '__main__':
    run()    

#!/usr/bin/env python

import os # This module provides a portable way of using operating system dependent functionality
EXAMPLE_DESCRIPTIVE_NAME = 'Example01'
EXAMPLE_AUTHOR = 'Pascual'
EXAMPLE_DIR = os.path.dirname(__file__) # Return the directory name of pathname path
GUIinclude = True

from datetime import date # The datetime module supplies classes for manipulating dates and times

# This module provides you easy access to your config/settings properties from all your python modules, it supports normal and lazy initialization for each property
from carputils import settings
from carputils import tools
from carputils import mesh
from carputils.carpio import txt

# Library for the Python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level 
# mathematical functions to operate on these arrays.
import numpy as np

def parser():
    parser = tools.standard_parser() #Generate a standard argument parser for collection of common options
    group  = parser.add_argument_group('experiment specific options') 
    group.add_argument('--duration',
                        type = float,
                        default = 2000.,
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
    """
    Generate name of top level output directory.
    """
    today = date.today()
    return '{}_basic_{}'.format(today.isoformat(), args.duration)

@tools.carpexample(parser, jobID)
def run(args, job):
   
    # Block which has a size of 5cm x 5cm x 1mm and a spatial resolution of 0.2 mm
    # Number of points (50/0.2)*(50/0.2)*(1/0.2) = 312500 points created that simulates the tissue
    #geom = mesh.Block(size=(50, 50, 1), resolution=0.2) #It is created an absolute path where it is stored a file of the mesh

    # Set fibre angle to 0, sheet angle to 0
    #geom.set_fibres(0, 0, 90, 90)

    # Generate and return base name
    #meshname = mesh.generate(geom)
    #meshname = '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/tissue_ints_10_10/tissue_ints_10_10'
    meshname = '/Volumes/yellow/Benutzer/pt732/tissue_ints_10_f1'

    # Mesh points to obtain the signal sources of the simulation. It must be taken into account that each source signal spans a region of the tissue (mesh created)
    # Number of points (50/1)*(50/1) = 2500 points created that simulated the signal sources
    points = mesh.Block(size=(45, 45, 0), resolution=3) #It is created an absolute path where it is stored a file of the mesh 
    points.set_fibres(0, 0, 0, 0)
    meshname2 = mesh.generate(points)
    
    # The geometry of the mesh to obtain the signal sources is created in (0,0,0) and it is wanted to move the position in the z axis on the surface of the mesh
    # called "geom" that simulates the tissue. This due to the fact that it is wanted to obtain the signals on the surface of the tissue
    # In order to do that, it will be added the position of the mesh of source points in the z axis to 0.55 mm = 550 µm, getting a correct position in the z axis
    dir_name = points.mesher_opts(meshname2)[25]
    points_text = np.loadtxt(dir_name + '.pts',skiprows = 1) #It is loaded the file where it is stored the geometry of the source points (.pts file)
    points_text[:,0] = points_text[:,0] + 25000
    points_text[:,1] = points_text[:,1] + 25000
    points_text[:,2] = points_text[:,2] + 2300 #It is added the position of the 3rd column of the file to 550 µm
    #print(points_text)

    #Save the file of the mesh points fixed 
    print(dir_name +'.pts')
    f=open(dir_name +'.pts','w')
    #print(str(points_text.shape[0]).join('\n'))
    f.write('{}\n'.format(points_text.shape[0])) #openCarp
    np.savetxt(f, points_text, delimiter=' ')
    f.close()

    #phie_rec_ptf = dir_name #"Defines the basename for the phie recovery file.
    #phie_recovery_file = electrograms #"Defines the file name used to output recovered extracellular potentials.
    #phie_rec_meth = 1 #"Defines the method for recovering phie with monodomain runs."

    egm = ['-phie_rec_ptf', dir_name ,
    '-phie_recovery_file', 'electrograms',
    '-phie_rec_meth', 1]


    # Query for element labels
    #_, etags, _ = txt.read(meshname + '.elem')
    #etags = np.unique(etags)

    #Number to identify different regions of the mesh (intracellular and extracellular model)
    #IntraTags = etags[etags != 0].tolist() #Returns a list of the values different equal to 0
    #ExtraTags = etags.tolist().copy() #Copy the list

    # Define the geometry of the stimulus at one end of the block
    # Units are um

    #stimtype = 0 (modelo intracelular)... diferentes modelos

    #Stimulus options
    stim = ['-num_stim',             2,
            '-stimulus[0].name',    'S1',
            '-stimulus[0].stimtype', 0, 
            '-stimulus[0].strength', args.S1_strength,
            '-stimulus[0].duration', args.S1_dur ,
            '-stimulus[1].name',    'S2',
            '-stimulus[1].stimtype', 0,
            '-stimulus[1].strength', args.S1_strength,
            '-stimulus[1].duration', args.S1_dur ,
            '-stimulus[0].x0', 0,
            '-stimulus[0].y0', 0,
            '-stimulus[0].z0', 800,
            '-stimulus[0].xd', 0,
            '-stimulus[0].yd', 50000,
            '-stimulus[0].zd', 2000,
            '-stimulus[1].start', 0,
            '-stimulus[1].x0', 0,
            '-stimulus[1].y0', 0,
            '-stimulus[1].z0', 800,
            '-stimulus[1].xd', 20000,
            '-stimulus[1].yd', 20000,
            '-stimulus[1].zd', 2000,
            '-stimulus[1].start', 192]

    #Generate a boundary condition definition

    #def carputils.mesh.block.block_boundary_condition(block, entity, index, coord, lower = True, bath = False, verbose = False)	
        #block	Block The geometry to generate boundary conditions for
        #entity	str The openCARP boundary condition type to use (e.g. 'stimulus')
        #index	int The boundary condition index
        #coord	str The direction perpendicular to which the BC is generated. One of ('x', 'y', 'z').
        #lower	bool, optional True if to generate BC at lower bound of coord (default), False for upper bound
        #bath	bool, optional Generate BC at bounds of tissue+bath if True, tissue only otherwise	
    #electrode = mesh.block_boundary_condition(geom, 'stimulus', 0, 'x', True)

    #Ionic setup
    imp = ['-num_imp_regions',3,
    '-imp_region[0].im','COURTEMANCHE',
    '-imp_region[0].num_IDs',1,
    '-imp_region[0].ID',1,
    '-imp_region[0].im_param', "g_CaL -55% , g_K1 +100% , blf_i_Kur -50% , g_to -65% , g_Ks +100% , maxI_pCa +50% , maxI_NaCa +60% , g_Kr *1.6 ",
    '-imp_region[1].im','COURTEMANCHE',
    '-imp_region[1].num_IDs',1,
    '-imp_region[1].ID',2,
    '-imp_region[1].im_param', "g_CaL -55% , g_K1 +100% , blf_i_Kur -50% , g_to -65% , g_Ks +100% , maxI_pCa +50% , maxI_NaCa +60% , g_Kr *1.6 ",
    '-imp_region[2].im','MacCannell_Fb',
    '-imp_region[2].num_IDs',1,
    '-imp_region[2].ID',3,
    ]
    
    #3 // COUNTERMANCHE remodelado //

    #Velocidad de conducción, región con ID = 0 se asignan velocidades de conducción en diferentes direcciones (longitudinal,transversal y normal)
    cv=['-num_gregions',3,
    '-gregion[0].num_IDs',1,
    '-gregion[0].ID',1,
    '-gregion[0].g_il',0.1438,
    '-gregion[0].g_it',0.1438,
    '-gregion[0].g_in',0.1438,
    '-gregion[1].num_IDs',1,
    '-gregion[1].ID',2,
    '-gregion[1].g_il',0.0719,
    '-gregion[1].g_it',0.0719,
    '-gregion[1].g_in',0.0719,
    '-gregion[2].num_IDs',1,
    '-gregion[2].ID',3,
    '-gregion[2].g_il',0.01438,
    '-gregion[2].g_it',0.01438,
    '-gregion[2].g_in',0.01438,
    ] #3 //

    #external/carputils/bin/tuneCV --model COURTEMANCHE --velocity 0.43 --sourceModel monodomain --converge True --dt 50
    #Conduction velocity: 0.4278 m/s [gi=0.1438, ge=0.5164, gm=0.1125]

    # Get basic command line, including solver options (ionic setup can be in other file)
    cmd = tools.carp_cmd() #Constructs a list of command line arguments. This will automatically 
                                                                 #include the loading of the correct options for the specified solvers

    cmd += ['-simID',    job.ID,
            '-meshname', meshname,
            '-dt',       50,
            '-tend',     args.duration]

    #cmd += ['--velocity', 43] #Velocity

    #cmd += tools.gen_physics_opts(ExtraTags=ExtraTags, IntraTags=IntraTags)
    cmd += stim + imp + cv + egm

    if args.visualize:
        cmd += ['-gridout_i', 3]
        cmd += ['-gridout_e', 3]
        cmd += ['-spacedt', 0.1]


    # Run simulation 
    job.carp(cmd) 

    
if __name__ == '__main__':
    run()    

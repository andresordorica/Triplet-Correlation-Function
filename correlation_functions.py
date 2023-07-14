
# Triplet-Correlation-Function
#    This function is based upon the methodology described by McNeil, Madden, Haymet and Rice:  
#    --The Journal of Chemical Physics 78, 388 (1983); doi: 10.1063/1.444514 
#    With further testing by Dhabal et al.
#    -- Phys. Chem. Chem. Phys., 2017,19, 3265-3278 
#    -- J. Chem. Phys. 7 November 2014; 141 (17): 174504. https://doi.org/10.1063/1.4898755
# The functions provided in this code employ the Kirkwood-Superposition-Principle for determining the triplet correlation:
#    g(r,s,t): Where r is the distance between p1 to p2, s the distance between p1 and p3, and t the distance between p2 and p3.
#    In this case p1 is the central particle with neighbors p1 and p2.
#    Both the self and non-self correlation functions assume r=s, with t being determined by the angle formed between the vectors p1:p2 and p1:p3
#    Consequently the distance "t" is determined based upon the cosine rule as t**2 = r**2 *( 2 *( 1 -cos(theta)))
#    g(r,s,t) whre r=s=t implies an equilateral configuration with a theta = 60 degrees.
#Authors: Andres Ordorica Fernandez
#This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

import numpy as np
import time
import sys
import os
import itertools
from utilities import *


sys.path.append(os.getcwd()+"/triplet_correlation_function")
from triplet_correlation_function import _tcf_self, _tcf_

'''
Fast_tcf_self and fast_tcf_ are ready to handle 

'''

def fast_tcf_self(u, atom_name_1,stride = 1000, theta = 60, plot=False, demo=True):
    '''
    Function to determine the triple correlation function for same atom types
    u: MDA analysis trajectory
    atom_name_1: (str) Selection for slicing trajectory, f.ex: "name O" (Selects all oxygen atoms from trajectory)
    stride: (int) number of frames to take ( 0 being every frame)
    theta: (int or float) angle in degrees to determine
    '''
    if demo ==False:
        #Example using a trajectory defined by the variable "u" loaded using MDAnalysis.
        # Select the target group, in this case the oxygen atom in water 
        atom_name_1 = "name O"
        atom_g1 = u.select_atoms('{}'.format(atom_name_1))
        traj1 =[]
        #Loop every nth frame of the trajectory extracting the target atoms positions
        #Convert from angstroms to nm
        #For computational efficiency an example block for selecting only 600 oxygen atoms
        for i, ts in enumerate(u.trajectory): 
            if i % stride ==0:
                coordinates_1 = np.array(atom_g1.positions) / 10
                # Check the number of rows in coordinates_1
                n_rows = coordinates_1.shape[0]
                # Select the first 500 rows if there are more than 500 rows
                if n_rows > 600:
                    coordinates_1 = coordinates_1[:600]
                traj1.append(coordinates_1)

        positions1 = np.array(traj1) # Array of positions with shape (M,N,3) where M = Number of frames, N = Number of atoms and 3 (x,y,z positions)
    
    else:
        #Performance and testing code
        N = 80  # Number of rows in each matrix ( In the trajectory data N = n_atoms of type a)
        M = 10  # Number of matrixes in the 3D array ( In trajectory data M = number of frames )
        #Create a 3D array of shape (M, N, 3) with random numbers
        positions1 = np.random.rand(M, N, 3)
        
    valid_distances, valid_angles = _tcf_self(positions1, tolerance = 0.010, angle = theta)
    if plot==True:
        plot_(valid_distances, theta, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), name_1 = atom_name_1.split()[1], name_2 = atom_name_1.split()[1])

    
    return valid_distances, valid_angles


def fast_tcf_(u, atom_name_1,atom_name_2,stride = 1000, theta = 60, plot=False, demo=True):
    if demo == False:
        '''
        Function to determine the triple correlation function for DIFFERENT atom types
        u: MDA analysis trajectory
        atom_name_1: (str) Selection for slicing trajectory, f.ex: "name Na" (Selects all sodium atoms from trajectory)
        atom_name_2: (str) Selection for slicing trajectory, f.ex: "name O" (Selects all oxygen atoms from trajectory)
        stride: (int) number of frames to take ( 0 being every frame)
        theta: (int or float) angle in degrees to determine
        '''
        
        # Select the groups where atom_g1 will be the central atoms 
        atom_name_1 = "name Na"
        atom_name_2 = "name O"
        atom_g1 = u.select_atoms('{}'.format(atom_name_1)) ############## ION
        atom_g2 = u.select_atoms('{}'.format(atom_name_2))
        traj1 =[]
        #Loop every nth frame of the trajectory extracting the target atoms positions
        #Convert from angstroms to nm
        #For computational efficiency an example block for selecting only 2000 oxygen atoms
        # In this block the positions of the central atom are stored (sodium) separately from the surrounding 
        #Atoms (Oxygen), after creating two arrays of position the arrays arestacked vertically to create a single matrix with 
        #both atoms in it.
        for i, ts in enumerate(u.trajectory):
            if i % stride ==0:
                coordinates_1 = np.array(atom_g1.positions) / 10 ################# ION
                if i ==0:
                    A = coordinates_1.shape[0]
                coordinates_2 = np.array(atom_g2.positions) / 10 ####################SOLVENT
                # Check the number of rows in coordinates_2
                n_rows2 = coordinates_2.shape[0]
                # Select the first 2000 rows if there are more than 2000 rows
                if n_rows2 > 2000:
                    coordinates_2 = coordinates_2[:2000]
                ###########################################################
                ### Stack the two arrays, where top (coordinates 1) is the ion and bottom (coordinates 2) is solvent
                combined_positions = np.concatenate((coordinates_1, coordinates_2), axis=0)
                if i ==0:
                    B = coordinates_2.shape[0]
                    print("Solvent array length: {}".format(str(A)))
                    C = combined_positions.shape[0]
                    print("total array length: {}".format(str(C)))
                    
                traj1.append(combined_positions)

        positions1 = np.array(traj1)
    
    else:
        n = 50
        a_ = 10
        M = 5
        positions = []
        for i in range(M):
            positions_atom_name_1 = np.random.rand(a_, 3)  # Example random positions data
            A = positions_atom_name_1.shape[1] # Take number of rows from target positions array (number of rows = number of atoms)
            positions_atom_name_2 = np.random.rand(n, 3)  # Example random positions data
            combined_positions = np.vstack((positions_atom_name_1, positions_atom_name_2))
            positions.append(combined_positions)
        positions1 = np.array(positions)
    
    valid_distances, valid_angles = _tcf_(positions1 , A, tolerance = 0.010, angle = theta)
    if plot==True:
        plot_(valid_distances, theta, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), name_1 = atom_name_1.split()[1], name_2 = atom_name_2.split()[1])

   
    return valid_distances, valid_angles

'''
In this blocks "u" should be replaced with an MDAnalysis trajectory
'''

#valid_distances_A_A, valid_angles_A_A  = fast_tcf_self("u", "name O", stride =1000, theta = 60, plot=False)


#valid_distances_A_B, valid_angles_A_B  = fast_tcf_("u", "name Na", "name O", stride =1000, theta = 60, plot=False)



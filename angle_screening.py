

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
from correlation_functions import *
from utilities import *


sys.path.append(os.getcwd()+"/triplet_correlation_function")
from triplet_correlation_function import _tcf_self, _tcf_



# Performing the analysis over a range of angles

def Angle_Screening(trajectory, 
                    atom_name_1 = "name A",
                    atom_name_2 =None, 
                    angle_range = [0, 180], 
                    angle_step = 20, 
                    stride = 10, 
                    plot_individual = False , 
                    demo_flag = True , 
                    chunk_size = 100, 
                    chunk_solvent_size = 100):
    angle_list = [i for i in range(angle_range[0], angle_range[1]+1 ,angle_step)]
    valid_distance_ =[]
    unbinned_ =[]
    var = demo_flag
    s = stride
    if atom_name_2 is None:
        for theta_ in angle_list:
            valid_distances_A_A, valid_angles_A_A  = fast_tcf_self(trajectory, atom_name_1, stride =s, theta = theta_, plot=plot_individual, demo = var, chunk = chunk_size)
            if demo_flag == False:
                unit_cell = (trajectory.trajectory.ts.triclinic_dimensions)/10
                unit_cell_volumes_ = (np.abs(np.linalg.det(unit_cell)))
                unit_cell_volumes = np.repeat(unit_cell_volumes_, int(trajectory.trajectory.n_frames))
                arr, r_ = hist(valid_distances_A_A, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), ucv = unit_cell_volumes)
            else:
                arr, r_ = hist(valid_distances_A_A, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), ucv = None)

            valid_distance_.append(arr)
            unbinned_.append(valid_distances_A_A)
    else:
        for theta_ in angle_list:
            valid_distances_A_B, valid_angles_A_B  = fast_tcf_(trajectory, atom_name_1,atom_name_2, stride =s, theta = theta_, plot=plot_individual, demo = var, chunk= chunk_size, chunk_solvent = chunk_solvent_size)
            if demo_flag == False:
                unit_cell = (trajectory.trajectory.ts.triclinic_dimensions)/10
                unit_cell_volumes_ = (np.abs(np.linalg.det(unit_cell)))
                unit_cell_volumes = np.repeat(unit_cell_volumes_, int(trajectory.trajectory.n_frames))
                arr, r_ = hist(valid_distances_A_B, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), ucv = unit_cell_volumes)
            else:
                arr, r_ = hist(valid_distances_A_B, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), ucv = None)
           
            valid_distance_.append(arr)
            unbinned_.append(valid_distances_A_B)


    plot_1 = lines_3D(valid_distance_, r_, angle_list, atom_name_1, atom_name_2)
    plot_2 = go_Mesh(valid_distance_, r_, angle_list, atom_name_1, atom_name_2)
    plot_3 = go_Histogram(np.array(unbinned_), r_, angle_list, atom_name_1, atom_name_2, angle_range = [0, 180], angle_step = 20, color_map = True, go_ = True)


    

     

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

def Angle_Screening(trajectory,atom_name_1 = "name A", atom_name_2 =None, angle_range = [0, 180], angle_step = 20, stride = 1000, plot_individual = False ):
    angle_list = [i for i in range(angle_range[0], angle_range[1]+1 ,angle_step)]
    valid_distance_ =[]
    unbinned_ =[]
    if atom_name_2 is None:
        for theta_ in angle_list:
            valid_distances_A_A, valid_angles_A_A  = fast_tcf_self(trajectory, atom_name_1, stride =stride, theta = theta_, plot=plot_individual)
            arr, r_ = hist(valid_distances_A_A, bin_width  = 0.01, r_range = np.array([0.2, 0.5]))
            valid_distance_.append(arr)
            unbinned_.append(valid_distances_A_A)
    else:
        for theta_ in angle_list:
            valid_distances_A_B, valid_angles_A_B  = fast_tcf_(trajectory, atom_name_1,atom_name_2, stride =stride, theta = theta_, plot=plot_individual)
            arr, r_ = hist(valid_distances_A_B, bin_width  = 0.01, r_range = np.array([0.2, 0.5]))
            valid_distance_.append(arr)
            unbinned_.append(valid_distances_A_B)


    plot_1 = lines_3D(valid_distance_, r_, angle_list, atom_name_1, atom_name_2)
    plot_2 = go_Mesh(valid_distance_, r_, angle_list, atom_name_1, atom_name_2)
    plot_3 = go_Histogram(np.array(unbinned_), r_, angle_list, atom_name_1, atom_name_2, angle_range = [0, 180], angle_step = 20, color_map = True, go_ = True)


    

     

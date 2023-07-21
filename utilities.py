
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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.io as pio
import warnings
warnings.filterwarnings("ignore")

sys.path.append(os.getcwd()+"/triplet_correlation_function")
from triplet_correlation_function import _tcf_self, _tcf_



def plot_(array, theta, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), name_1 = "None", name_2 = "None", ucv = None):
        fig, ax = plt.subplots()
        n_bins = int((r_range[1] - r_range[0]) / bin_width)
        g_rst, edges = np.histogram(array, range=r_range, bins=n_bins)
        r = 0.5 * (edges[1:] + edges[:-1]) 
        V = (4/3)* (np.pi) * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        if ucv is None:
             norm = 3 *V 
             print("Using demo, no trajectory unit cell volumes for normalization")
        else:
            norm =  3 * np.sum(1.0 / ucv) * V
       
        try:
            g_rst = g_rst.astype(np.float64) / norm  # From int64.  
        except:
            g_rst = g_rst/ norm
        

       
        ax.set_xlabel(r'$r$ (nm)')
        if theta == 60:
            if name_2 is not None:
                ax.set_ylabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,r)'.format(name_2, name_1, name_2))
            else:
                ax.set_ylabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,r)'.format(name_1, name_1, name_1))
                
        elif theta != 60:
            if name_2 is not None:
                ax.set_ylabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(name_2, name_1, name_2))
            else:
                ax.set_ylabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(name_1, name_1, name_1))
        ax.set_title("Theta: {} degrees".format(str(theta)))
        ax.plot(r, g_rst, '-')
        plt.show()

def hist(array, bin_width  = 0.01, r_range = np.array([0.2, 0.5]), ucv = None):
    n_bins = int((r_range[1] - r_range[0]) / bin_width)
    g_rst, edges = np.histogram(array, range=r_range, bins=n_bins)
    r = 0.5 * (edges[1:] + edges[:-1]) 
    '''
        #load md traj tracjectory to extract the unitcell volumes and normalize the 
        #triplet correltion function
        V = (4/3)* (np.pi) * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        norm =  3 * np.sum(1.0 / trj.unitcell_volumes) * V
    '''
    V = (4/3)* (np.pi) * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    if ucv is None:
        norm = 3 *V 
        print("Using demo, no trajectory unit cell volumes for normalization")
    else:
        norm =  3 * np.sum(1.0 / ucv) * V


    try:
        g_rst = g_rst.astype(np.float64) / norm  # From int64.  
    except:
        g_rst = g_rst/ norm
    
    return g_rst, r
        
def lines_3D(binned_distances, r_, angle_list, atom_name_1, atom_name_2):
    # Create the figure and axes
    fig = plt.figure()
    axs = fig.add_subplot(111, projection='3d')
    matrix = np.empty((0, 3))
    # Iterate through the lines
    for i, dis in enumerate(binned_distances):
        # Generate data for each line
        t_T = (np.array(r_)**2) * 2 * (1 - np.cos(np.radians(angle_list[i])))
        t = [angle_list[i]]
        # Repeat the values in t to match the desired length
        t_ = np.tile(t, int(np.ceil(len(t_T) / len(t))))[:len(t_T)]

        #
        x = r_
        y = dis
        z = t_
        # Create a temporary matrix for the current iteration
        temp_matrix = np.concatenate((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]), axis=1)
        # Append the temporary matrix to the main matrix
        matrix = np.vstack((matrix, temp_matrix))
        axs.plot(x,z,y)
        axs.view_init(27, 21)  # Rotate 90 degrees around x-axis, -90 degrees around y-axis
    
   
    axs.set_xlabel(r'$r$ (nm)')
    axs.set_ylabel(r'$\theta$ (degrees) ')   
    if atom_name_2 is not None:
        axs.set_zlabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(atom_name_2.split()[-1], atom_name_1.split()[-1], atom_name_2.split()[-1]))
    else:
        axs.set_zlabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(atom_name_1.split()[-1], atom_name_1.split()[-1], atom_name_1.split()[-1]))
    plt.show()
       
def go_Mesh(binned_distances, r_, angle_list, atom_name_1, atom_name_2, color_map = True):
    matrix = np.empty((0, 3))
    for i, dis in enumerate(binned_distances):
        # Generate data for each line
        t_T = (np.array(r_)**2) * 2 * (1 - np.cos(np.radians(angle_list[i])))
        t = [angle_list[i]]
        # Repeat the values in t to match the desired length
        t_ = np.tile(t, int(np.ceil(len(t_T) / len(t))))[:len(t_T)]
        x = r_
        y = dis
        z = t_
        # Create a temporary matrix for the current iteration
        temp_matrix = np.concatenate((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]), axis=1)
        # Append the temporary matrix to the main matrix
        matrix = np.vstack((matrix, temp_matrix))

    # Data for three-dimensional scattered points
    x = matrix[:,0]
    y = matrix[:,1]
    z = matrix[:,2]
    # Normalize the y values to the range [0, 1]
    y_normalized = (y - y.min()) / (y.max() - y.min())
    fig = go.Figure(data=[go.Mesh3d(x=x, y=z, z=y, intensity=y_normalized, opacity=0.40, colorscale='Viridis')])
    # Update the layout with the axis labels
    if atom_name_2 is not None:
        fig.update_layout(
        scene=dict(
            xaxis=dict(title=r'r (nm)'),
            yaxis=dict(title=r'Theta (degrees)'),
            zaxis=dict(title=r'g_{}-{}-{} (r,r,s)'.format(atom_name_2.split()[-1],atom_name_1.split()[-1],atom_name_2.split()[-1]))
            )
            )
    else:
        fig.update_layout(
        scene=dict(
            xaxis=dict(title=r'r (nm)'),
            yaxis=dict(title=r'Theta (degrees)'),
            zaxis=dict(title=r'g_{}-{}-{} (r,r,s)'.format(atom_name_1.split()[-1],atom_name_1.split()[-1],atom_name_1.split()[-1]))
            )
            )

    fig.show()

def hist_3d(matrix,bin_width  = 0.01, r_range = np.array([0.2, 0.5]), angle_range = [0, 180], angle_step = 40):
    # Extract the x, y, and z columns from the matrix
    x = matrix[:, 0]
    z = matrix[:, 1]
    #X bins
    bins_x =np.linspace(r_range[0], r_range[1], int((r_range[1] - r_range[0]) / bin_width))
    bins_z =np.linspace(angle_range[0], angle_range[1], int((angle_range[1] - angle_range[0]) / angle_step))
    # Compute the 2D histogram
    hist, x_edges, y_edges = np.histogram2d(matrix[:, 0], matrix[:, 1], bins=[bins_x, bins_z])
    # Create a meshgrid for the bin edges
    X_edges, Y_edges = np.meshgrid(x_edges[:-1], y_edges[:-1])

    return hist, X_edges, Y_edges
    
def go_Histogram(binned_distances, r_, angle_list, atom_name_1, atom_name_2,angle_range = [0, 180], angle_step = 20, color_map = True, go_ = True):
    matrix = np.empty((0, 2))
    for i, dis in enumerate(binned_distances):
        # Generate data for each line
        t = [angle_list[i]]
        # Repeat the values in t to match the desired length
        t_ = np.tile(t, int(np.ceil(len(dis) / len(t))))[:len(dis)]
        y = dis
        z = t_
        # Create a temporary matrix for the current iteration
        temp_matrix = np.concatenate((y[:, np.newaxis], z[:, np.newaxis]), axis=1)
        # Append the temporary matrix to the main matrix
        matrix = np.vstack((matrix, temp_matrix))


    hist, X_edges, Y_edges = hist_3d(matrix,bin_width  = 0.01, r_range = np.array([0.2, 0.5]), angle_range = [0, 180], angle_step = 40)
    # Create a figure and axes
    if go_ ==False:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Plot the 3D histogram
        ax.plot_surface(X_edges, Y_edges, hist.T, cmap='viridis')
        # Set labels and title
        ax.set_xlabel(r'$r$ (nm)')
        ax.set_ylabel(r'$\theta$ (degrees) ')   
        if atom_name_2 is not None:
            ax.set_zlabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(atom_name_2.split()[-1], atom_name_1.split()[-1], atom_name_2.split()[-1]))
        else:
            ax.set_zlabel(r'$g_{{\mathrm{{{}-{}-{}}}}}$ (r,r,s)'.format(atom_name_1.split()[-1], atom_name_1.split()[-1], atom_name_1.split()[-1]))

        plt.title('3D plot of 2D histogram')
        plt.show()
    if go_ == True:
        fig = plt.figure()
        # Create the surface plot using Plotly
        fig = go.Figure(data=[go.Surface(x=X_edges, y=Y_edges, z=hist.T,  colorscale='Viridis', opacity=0.4)])

        # Set the axis labels and title
        fig.update_layout(scene=dict(
                            xaxis_title='r (nm)',
                            yaxis_title=' Theta (degrees)',
                            zaxis_title='Frequency'),
                        title='3D plot of 2D histogram')

        # Show the plot
        fig.show()





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

from utilities import * 

# Performing the analysis over a range of angles    

    
plot_3d = Angle_Screening("u",atom_name_1 = "name O", atom_name_2 =None, angle_range = [0, 180], angle_step = 20, stride = 1000 )
    
          
     

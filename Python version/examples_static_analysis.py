"""
-------------------------------------------------------------------------------
FEM Bar Elements Problems: Static Analysis
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt 
from FEM_matrices import FEM_matrices
from static_analysis import static_analysis
from plot_structure import plot_structure


# Example 1: Plane frame ------------------------------------------------------
# Input data
nodes = pd.read_csv('Input_planeframe_nodes.txt',delim_whitespace=True, header = None,dtype = 'float')
nodes.fillna(0, inplace=True)  
nodes = nodes.values 

bars = pd.read_csv('Input_planeframe_bars.txt',delim_whitespace=True, header = None,dtype = 'float')
bars.fillna(0, inplace=True)   
bars = bars.values 

element_type = 'plane frame'

# Figure
plt.figure(1)
plot_structure(nodes,bars,element_type)

# FEM Matrices  
K_pf, _, Kr_pf, _ = FEM_matrices(nodes,bars,element_type)    

# Solution: Nodal Displacements and Forces
disp_pf, forces_pf = static_analysis(nodes,element_type,K_pf,Kr_pf)

# Example 2: Truss structure --------------------------------------------------
# Input data
nodes = pd.read_csv('Input_truss_nodes.txt',delim_whitespace=True, header = None,dtype = 'float')
nodes.fillna(0, inplace=True)  
nodes = nodes.values 

bars = pd.read_csv('Input_truss_bars.txt',delim_whitespace=True, header = None,dtype = 'float')
bars.fillna(0, inplace=True)   
bars = bars.values 

element_type = 'plane truss'

# Figure
plt.figure(2)
plot_structure(nodes,bars,element_type)

# FEM Matrices  
K_truss, _, Kr_truss, _ = FEM_matrices(nodes,bars,element_type)     

# Solution: Nodal Displacements and Forces
disp_truss, forces_truss = static_analysis(nodes,element_type,K_truss,Kr_truss)
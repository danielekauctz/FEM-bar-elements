"""
-------------------------------------------------------------------------------
FEM Bar Elements Problems: Modal Parameters
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt 
from FEM_matrices import FEM_matrices
from modal_analysis import modal_analysis
from plot_structure import plot_structure

# Example 1: Timoshenko beam --------------------------------------------------
# Input data
nodes = pd.read_csv('Input_beam_nodes.txt',delim_whitespace=True, header = None,dtype = 'float')
nodes.fillna(0, inplace=True)  
nodes = nodes.values 

bars = pd.read_csv('Input_beam_bars.txt',delim_whitespace=True, header = None,dtype = 'float')
bars.fillna(0, inplace=True)   
bars = bars.values 

element_type = 'timoshenko beam'

# Figure
plt.figure(1)
plot_structure(nodes,bars,element_type)

# FEM Matrices  
_, _, Kr_beam, Mr_beam = FEM_matrices(nodes,bars,element_type)    

# Modal parameters 
fn_beam, _, phi_beam = modal_analysis(Kr_beam, Mr_beam)


# Example 2: Plane frame ------------------------------------------------------
# Input data
nodes = pd.read_csv('Input_planeframe_nodes.txt',delim_whitespace=True, header = None,dtype = 'float')
nodes.fillna(0, inplace=True)  
nodes = nodes.values 

bars = pd.read_csv('Input_planeframe_bars.txt',delim_whitespace=True, header = None,dtype = 'float')
bars.fillna(0, inplace=True)   
bars = bars.values 

element_type = 'plane frame'

# Figure
plt.figure(2)
plot_structure(nodes,bars,element_type)

# FEM Matrices  
_, _, Kr_pf, Mr_pf = FEM_matrices(nodes,bars,element_type)    

# Modal parameters 
fn_pf, _, phi_pf = modal_analysis(Kr_pf, Mr_pf)

# Example 3: Truss structure --------------------------------------------------
# Input data
nodes = pd.read_csv('Input_truss_nodes.txt',delim_whitespace=True, header = None,dtype = 'float')
nodes.fillna(0, inplace=True)  
nodes = nodes.values 

bars = pd.read_csv('Input_truss_bars.txt',delim_whitespace=True, header = None,dtype = 'float')
bars.fillna(0, inplace=True)   
bars = bars.values 

element_type = 'plane truss'

# Figure
plt.figure(3)
plot_structure(nodes,bars,element_type)

# FEM Matrices  
_, _, Kr_truss, Mr_truss = FEM_matrices(nodes,bars,element_type)    

# Modal parameters 
fn_truss, _, phi_truss = modal_analysis(Kr_truss, Mr_truss)
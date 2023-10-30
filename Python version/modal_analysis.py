"""
-------------------------------------------------------------------------------
NATURAL FREQUENCIES AND MODE SHAPES
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
Parameters:
 Kr, Mr: stiffness and mass matrices with boundary conditions
 fn: natural frequencies (Hz)
 wn: natural frequencies (rad/s)
 phi: mode shapes
-------------------------------------------------------------------------------
"""
# Importing modules
import scipy.linalg as sc
import numpy as np

# Function
def modal_analysis(Kr, Mr):

    w2, phi = sc.eig(Kr, Mr)

    iw  = w2.argsort()
    w2  = w2[iw]
    phi = phi[:,iw]

    wn  = np.sqrt(np.real(w2)) 
    fn  = wn/(2*np.pi)
    
    for kk in range(len(phi)):
       phi[:,kk] = (abs(phi[:,kk])/max(abs(phi[:,kk])))*np.sign(np.real(phi[:,kk]));

    return fn, wn, phi
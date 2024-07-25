"""
-------------------------------------------------------------------------------
Static Analysis: Nodal Displacements and Forces
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
Input data:
 nodes: text file with information about structure nodes
      column 1 and 2: node X and Y coordinate
      column 3, 4 and 5: node boundary conditions (horizontal displacement, 
             vertical displacement, rotation (if = 1 movement is restricted)
      column 6, 7 and 8: acting forces in the node (horizontal force, vertical 
             force, bending moment
      obs., truss elements have no rotation or bending moment, so there
              is no columns 5 and 8
 bars: text file with information about structure bars
     column 1 and 2: N1, N2 = element nodes
     column 3: A = cross sectional area
     column 4: E = Elastic modulus
     column 5: I = moment of inertia
     column 6: ro = specific mass
 element_type: bar element type
                 'plane frame'
                 'timoshenko beam'
                 'plane truss'
 K: stiffness matrix
 Kr: stiffness matrix with boundary conditions
-------------------------------------------------------------------------------
Output data:
 disp: nodal displacements vector
 forces: nodal forces vector
-------------------------------------------------------------------------------
"""
import numpy as np

def static_analysis(nodes,element_type,K,Kr):
    
    if element_type == 'plane frame' or element_type == 'timoshenko beam':
        dof = 3
        maxdof = dof*np.size(nodes,0)
        F = np.zeros((maxdof,1))
        restr = np.zeros(dof*np.size(nodes,0))
        
        for i in range(np.size(nodes,0)):
            RX = nodes[i,2]
            RY = nodes[i,3]
            RZ = nodes[i,4]
            Fx = nodes[i,5]
            Fy = nodes[i,6]
            Mz = nodes[i,7]
            
            gl1 = int(dof*(i+1)-2)
            gl2 = int(dof*(i+1)-1)
            gl3 = int(dof*(i+1))
            
            if RX == 1:
                restr[gl1-1] = 1
            if RY == 1:
                restr[gl2-1] = 1
            if RZ == 1:
                restr[gl3-1] = 1
            if Fx != 0:
                F[gl1-1] = Fx
            if Fy != 0:
                F[gl2-1] = Fy
            if Mz != 0:
                F[gl3-1] = Mz
                
    elif element_type == 'plane truss':
        dof = 2
        maxdof = dof*np.size(nodes,0)
        F = np.zeros((maxdof,1))
        restr = np.zeros(dof*np.size(nodes,0))
        
        for i in range(np.size(nodes,0)):
            RX = nodes[i,2]
            RY = nodes[i,3]
            Fx = nodes[i,4]
            Fy = nodes[i,5]
            
            gl1 = int(dof*(i+1)-1)
            gl2 = int(dof*(i+1))
            
            if RX == 1:
                restr[gl1-1] = 1
            if RY == 1:
                restr[gl2-1] = 1
            if Fx != 0:
                F[gl1-1] = Fx
            if Fy != 0:
                F[gl2-1] = Fy
                
    ccnt = restr==1
    ccnt = ccnt.nonzero()
    ccnt = np.asarray(ccnt)
    Fr = np.delete(F,ccnt, axis = 0)
    
    disp_r = np.dot(np.linalg.inv(Kr), Fr)
    disp = np.zeros((maxdof,1))
    ccnt = restr==0
    ccnt = ccnt.nonzero()
    disp[ccnt,0] = disp_r.T
    
    forces = K@disp
            
    return disp, forces

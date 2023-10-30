"""
-------------------------------------------------------------------------------
STIFFNESS AND MASS MATRICES FOR THE FINITE ELEMENT METHOD (FEM)
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
Input data:
 nodes: text file with information about structure nodes
      column 1 and 2: node X and Y coordinate
      column 3, 4 and 5: node boundary conditions (horizontal displacement, 
             vertical displacement, rotation (if = 1 movement is restricted)
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
-------------------------------------------------------------------------------
Output data:
 K: stiffness matrix
 M: consistent mass matrix
 Kr, Mr: stiffness and mass matrices with boundary conditions
-------------------------------------------------------------------------------
"""
# Importing modules
import numpy as np

# Function
def FEM_matrices(nodes,bars,element_type):
    
    if element_type == 'timoshenko beam' or element_type == 'plane frame':
        dof = 3  # degree of freedom per node
    elif element_type == 'plane truss':
        dof = 2
        
    # Elements properties
    L = np.zeros(np.size(bars,0))
    seno = np.zeros(np.size(bars,0))
    coss = np.zeros(np.size(bars,0))
    for i in range(np.size(bars,0)):
        N1 = int(bars[i,0])
        N2 = int(bars[i,1])

        x1 = nodes[N1-1,0]
        y1 = nodes[N1-1,1]
        x2 = nodes[N2-1,0]
        y2 = nodes[N2-1,1]

        Lx = x2 - x1
        Ly = y2 - y1

        Ls = np.sqrt(Lx**2 + Ly**2)
        L[i] = Ls
        seno[i]  = Ly/Ls
        coss[i]  = Lx/Ls
        
    # Matrices
    maxdof = dof*np.size(nodes,0)

    K = np.zeros((maxdof,maxdof))
    M = np.zeros((maxdof,maxdof))

    for i in range(np.size(bars,0)):
        Ls = L[i]
        sen = seno[i]
        cos = coss[i]
        N1 = bars[i,0]
        N2 = bars[i,1]
        A = bars[i,2]
        E = bars[i,3]
        I = bars[i,4]
        ro = bars[i,5]
        
        if element_type == 'plane frame':
            K11 = (E*A)/Ls
            K22 = (12*E*I)/(Ls**3)
            K32 = (6*E*I)/(Ls**2)
            K23 = K32
            K33 = (4*E*I)/Ls
            K36 = (2*E*I)/Ls
            
            K1 = np.array([[ K11,    0,     0, -K11,    0,     0],
                            [    0,  K22,   K23,    0, -K22,   K32],  
                            [    0,  K32,   K33,    0, -K32,   K36], 
                            [ -K11,    0,     0,  K11,    0,     0],
                            [    0, -K22,  -K23,    0,  K22,  -K32],
                            [    0,  K32,   K36,    0, -K32,   K33]])
                          
            M1 = ((ro*A*Ls)/420)*np.array([[140,      0,       0,  70,      0,       0],
                                            [  0,    156,   22*Ls,   0,     54,  -13*Ls],  
                                            [  0,  22*Ls,  4*Ls**2,   0,  13*Ls, -3*Ls**2], 
                                            [ 70,      0,       0, 140,      0,       0],
                                            [  0,     54,   13*Ls,   0,    156,  -22*Ls],
                                            [  0, -13*Ls, -3*Ls**2,   0, -22*Ls,  4*Ls**2]]) 
            
            Rot = np.array([[ cos,  sen,    0,    0,   0,   0],          # Rotation matrix
                            [-sen,  cos,    0,    0,   0,   0],
                            [   0,    0,    1,    0,   0,   0],
                            [   0,    0,    0,  cos, sen,   0],
                            [   0,    0,    0, -sen, cos,   0],
                            [   0,    0,    0,    0,   0,   1]])
            
            Klr = np.dot(np.dot(Rot.T, K1), Rot)
            Mrr = np.dot(np.dot(Rot.T, M1), Rot)
            
            gl1 = int(dof*N1-2)
            gl2 = int(dof*N1-1)
            gl3 = int(dof*N1)
            gl4 = int(dof*N2-2)
            gl5 = int(dof*N2-1)
            gl6 = int(dof*N2)
            
            K[gl1-1:gl3, gl1-1:gl3] += Klr[0:3, 0:3]
            K[gl4-1:gl6, gl1-1:gl3] += Klr[3:6, 0:3]
            K[gl1-1:gl3, gl4-1:gl6] += Klr[0:3, 3:6]
            K[gl4-1:gl6, gl4-1:gl6] += Klr[3:6, 3:6]
            
            M[gl1-1:gl3, gl1-1:gl3] += Mrr[0:3, 0:3]
            M[gl4-1:gl6, gl1-1:gl3] += Mrr[3:6, 0:3]
            M[gl1-1:gl3, gl4-1:gl6] += Mrr[0:3, 3:6]
            M[gl4-1:gl6, gl4-1:gl6] += Mrr[3:6, 3:6]
            
        elif element_type == 'timoshenko beam':
            nu = 0.3                                 # Poisson’s coefficient
            k = 0.5                                  # Timoshenko’s shear factor               
            G = E/2/(1+nu)                           # Elastic shear modulus
            
            fi = (12/Ls**2)*(E*I/(k*G*A))
            C = E*I/((1+fi)*Ls**3)
            
            K1 = C*np.array([[(A/I)*(1+fi),    0,            0, -(A/I)*(1+fi),     0,            0],
                             [           0,   12,         6*Ls,             0,   -12,         6*Ls],  
                             [           0, 6*Ls, (4+fi)*Ls**2,             0, -6*Ls, (2-fi)*Ls**2],
                             [-(A/I)*(1+fi),   0,            0,  (A/I)*(1+fi),     0,            0],
                             [           0,  -12,        -6*Ls,             0,    12,        -6*Ls],
                             [           0, 6*Ls, (2-fi)*Ls**2,             0, -6*Ls, (4+fi)*Ls**2]])
                               
            M1 = ((ro*A*Ls)/(210*(1+fi)**2))*np.array([[70*(1+fi)**2,                            0,                             0, 35*(1+fi)**2,                            0,                             0],
                                                       [           0,       (70*fi**2)+(147*fi)+78,  ((35*fi**2)+(77*fi)+44)*Ls/4,            0,        (35*fi**2)+(63*fi)+27, -((35*fi**2)+(63*fi)+26)*Ls/4],  
                                                       [           0, ((35*fi**2)+(77*fi)+44)*Ls/4, ((7*fi**2)+(14*fi)+8)*Ls**2/4,            0, ((35*fi**2)+(63*fi)+26)*Ls/4,-((7*fi**2)+(14*fi)+6)*Ls**2/4], 
                                                       [35*(1+fi)**2,                            0,                             0, 70*(1+fi)**2,                            0,                             0],
                                                       [           0,        (35*fi**2)+(63*fi)+27,  ((35*fi**2)+(63*fi)+26)*Ls/4,            0,       (70*fi**2)+(147*fi)+78, -((35*fi**2)+(77*fi)+44)*Ls/4],
                                                       [           0,-((35*fi**2)+(63*fi)+26)*Ls/4,-((7*fi**2)+(14*fi)+6)*Ls**2/4,            0,-((35*fi**2)+(77*fi)+44)*Ls/4,((7*fi**2)+(14*fi)+8)*Ls**2/4]]) 
               
            M2 = ((ro*I)/(30*(1+fi)**2*Ls))*np.array([[0,               0,                           0, 0,              0,                           0],
                                                      [0,              36,             -((15*fi)-3)*Ls, 0,            -36,             -((15*fi)-3)*Ls],
                                                      [0, -((15*fi)-3)*Ls, ((10*fi**2)+(5*fi)+4)*Ls**2, 0, ((15*fi)-3)*Ls,  ((5*fi**2)-(5*fi)-1)*Ls**2],
                                                      [0,               0,                           0, 0,              0,                           0],
                                                      [0,             -36,              ((15*fi)-3)*Ls, 0,             36,              ((15*fi)-3)*Ls],
                                                      [0, -((15*fi)-3)*Ls,  ((5*fi**2)-(5*fi)-1)*Ls**2, 0, ((15*fi)-3)*Ls, ((10*fi**2)+(5*fi)+4)*Ls**2]])

            Mt = M1 + M2
            
            Rot = np.array([[ cos,  sen,    0,    0,   0,   0],
                            [-sen,  cos,    0,    0,   0,   0],
                            [   0,    0,    1,    0,   0,   0],
                            [   0,    0,    0,  cos, sen,   0],
                            [   0,    0,    0, -sen, cos,   0],
                            [   0,    0,    0,    0,   0,   1]])
            
            Klr = np.dot(np.dot(Rot.T, K1), Rot)
            Mrr = np.dot(np.dot(Rot.T, Mt), Rot)
            
            gl1 = int(dof*N1-2)
            gl2 = int(dof*N1-1)
            gl3 = int(dof*N1)
            gl4 = int(dof*N2-2)
            gl5 = int(dof*N2-1)
            gl6 = int(dof*N2)
            
            K[gl1-1:gl3, gl1-1:gl3] += Klr[0:3, 0:3]
            K[gl4-1:gl6, gl1-1:gl3] += Klr[3:6, 0:3]
            K[gl1-1:gl3, gl4-1:gl6] += Klr[0:3, 3:6]
            K[gl4-1:gl6, gl4-1:gl6] += Klr[3:6, 3:6]
            
            M[gl1-1:gl3, gl1-1:gl3] += Mrr[0:3, 0:3]
            M[gl4-1:gl6, gl1-1:gl3] += Mrr[3:6, 0:3]
            M[gl1-1:gl3, gl4-1:gl6] += Mrr[0:3, 3:6]
            M[gl4-1:gl6, gl4-1:gl6] += Mrr[3:6, 3:6]
            
        elif element_type == 'plane truss':
            Ke = (E*A)/Ls
            K1 = Ke*np.array([[ 1, 0, -1, 0],
                              [ 0, 0,  0, 0],  
                              [-1, 0,  1, 0], 
                              [ 0, 0,  0, 0]])
                      
            Me = (ro*A*Ls)/6
            M1 = Me*np.array([[2, 0, 1, 0],
                              [0, 2, 0, 1],
                              [1, 0,  2, 0], 
                              [0, 1,  0, 2]])

            Rot = np.array([[ cos, sen,   0,   0],
                            [-sen, cos,   0,   0],  
                            [   0,   0, cos, sen], 
                            [   0,   0,-sen, cos]])
            
            Klr = np.dot(np.dot(Rot.T, K1), Rot)
            Mrr = np.dot(np.dot(Rot.T, M1), Rot)
            
            gl1 = int(dof*N1-1)
            gl2 = int(dof*N1)
            gl3 = int(dof*N2-1)
            gl4 = int(dof*N2)

            K[gl1-1:gl2, gl1-1:gl2] += Klr[0:2, 0:2]
            K[gl3-1:gl4, gl1-1:gl2] += Klr[2:4, 0:2]
            K[gl1-1:gl2, gl3-1:gl4] += Klr[0:2, 2:4]
            K[gl3-1:gl4, gl3-1:gl4] += Klr[2:4, 2:4]
            
            M[gl1-1:gl2, gl1-1:gl2] += Mrr[0:2, 0:2]
            M[gl3-1:gl4, gl1-1:gl2] += Mrr[2:4, 0:2]
            M[gl1-1:gl2, gl3-1:gl4] += Mrr[0:2, 2:4]
            M[gl3-1:gl4, gl3-1:gl4] += Mrr[2:4, 2:4]
            
    # Boundary conditions
    restr = np.zeros(dof*np.size(nodes,0))
    for i in range(np.size(nodes,0)):
        if element_type == 'timoshenko beam' or element_type == 'plane frame':
            RX = nodes[i,2]
            RY = nodes[i,3]
            RZ = nodes[i,4]
            
            gl1 = int(dof*(i+1)-2)
            gl2 = int(dof*(i+1)-1)
            gl3 = int(dof*(i+1))
            
            if RX == 1:
                restr[gl1-1] = 1
            if RY == 1:
                restr[gl2-1] = 1
            if RZ == 1:
                restr[gl3-1] = 1
            
        elif element_type == 'plane truss':
            RX = nodes[i,2]
            RY = nodes[i,3]
            
            gl1 = int(dof*(i+1)-1)
            gl2 = int(dof*(i+1))
            
            if RX == 1:
                restr[gl1-1] = 1
            if RY == 1:
                restr[gl2-1] = 1
                
    ccnt = restr==1
    ccnt = ccnt.nonzero()
    Kr = np.delete(K,ccnt, axis = 0)  
    Kr = np.delete(Kr,ccnt, axis = 1)  
    Mr = np.delete(M,ccnt, axis = 0)  
    Mr = np.delete(Mr,ccnt, axis = 1)     
    
    return K, M, Kr, Mr

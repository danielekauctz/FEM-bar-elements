"""
-------------------------------------------------------------------------------
STRUCTURE FIGURE
Daniele Kauctz Monteiro (2023)
danielekauctz@hotmail.com
-------------------------------------------------------------------------------
Input data:
 nodes: text file with information about structure nodes
           column 1 and 2: node X and Y coordinate
           column 3, 4 and 5: node boundary conditions (horizontal displacement, 
             vertical displacement, rotation (if = 1 movement is restricted)
             obs., truss elements have no rotation, so there is no columns 
             5 and 8
 bars: text file with information about structure bars
           column 1 and 2: N1, N2 = element nodes         
 element_type: bar element type
                 'plane frame'
                 'timoshenko beam'
                 'plane truss'
-------------------------------------------------------------------------------
"""
import numpy as np
import matplotlib.pyplot as plt 

def plot_structure(nodes,bars,element_type):
    for i in range(np.size(bars,0)):
        N1 = int(bars[i,0])
        N2 = int(bars[i,1])

        x1 = nodes[N1-1,0]
        y1 = nodes[N1-1,1]
        x2 = nodes[N2-1,0]
        y2 = nodes[N2-1,1]
        
        x  = np.array([x1, x2])
        y  = np.array([y1, y2])
        
        xmax = max(nodes[:,0])
        if xmax == 0:
            xmax = 1
            
        ymax = max(nodes[:,1]) 
        if ymax == 0:
            ymax = 1
            
        plt.plot(x, y, c='black')
        plt.scatter(x, y, c='black', marker='o')
        plt.xlim([-0.2*xmax, xmax+(0.2*xmax)])
        plt.ylim([-0.1*ymax, ymax+(0.1*ymax)])
        
    for i in range(np.size(nodes,0)):
        X = nodes[i,0]
        Y = nodes[i,1]
        RX = nodes[i,2]
        RY = nodes[i,3]
            
        if RX == 1:                   # zero X-displacement boundary conditions
            plt.scatter(X, Y, 300, c='r', marker=5, zorder = -2)
            
        if RY == 1:                   # zero Y-displacement boundary conditions
            plt.scatter(X, Y, 300, c='r', marker=6, zorder = -2)
            
        if element_type == 'plane frame' or element_type == 'timoshenko beam':
            RZ = nodes[i,4]
            if RZ == 1:               # zero Z-displacement boundary conditions
                plt.scatter(X, Y, 800, c='r', marker='s', zorder = -2)
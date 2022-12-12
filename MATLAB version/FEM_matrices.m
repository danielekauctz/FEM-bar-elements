% -------------------------------------------------------------------------
% STIFFNESS AND MASS MATRICES FOR THE FINITE ELEMENT METHOD (FEM)
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com
% -------------------------------------------------------------------------
% Input data:
% nodes: text file with information about structure nodes
    ... column 1 and 2: node X and Y coordinate
    ... column 3, 4 and 5: node boundary conditions (horizontal displacement, 
%             vertical displacement, rotation (if = 1 movement is restricted)
% bars: text file with information about structure bars
   ... column 1 and 2: N1, N2 = element nodes
   ... column 3: A = cross sectional area
   ... column 4: E = Elastic modulus
   ... column 5: I = moment of inertia
   ... column 6: ro = specific mass
% element_type: bar element type
               ... 'plane frame'
               ... 'timoshenko beam'
               ... 'plane truss'
% -------------------------------------------------------------------------
% Output data:
% K: stiffness matrix
% M: consistent mass matrix
% Kr, Mr: stiffness and mass matrices with boundary conditions
% -------------------------------------------------------------------------

function [K,M,Kr,Mr] = FEM_matrices(nodes,bars,element_type)
    
    if strcmp(element_type,'timoshenko beam') || strcmp(element_type,'plane frame')
        dof = 3;  % degree of freedom per node
    elseif strcmp(element_type,'plane truss')
        dof = 2;
    end

    % Elements properties
    for i = 1:size(bars,1)
        N1 = bars(i,1);
        N2 = bars(i,2);

        x1 = nodes(N1,1);
        y1 = nodes(N1,2);
        x2 = nodes(N2,1);
        y2 = nodes(N2,2);

        Lx = x2 - x1;
        Ly = y2 - y1;

        Ls = sqrt(Lx^2 + Ly^2);
        L(i) = Ls;
        seno(i)  = Ly/Ls;
        coss(i)  = Lx/Ls;
    end

    % Matrices
    maxdof = dof*size(nodes,1);

    K = zeros(maxdof,maxdof);
    M = zeros(maxdof,maxdof);

    for i = 1:size(bars,1)
        Ls = L(i);
        sen = seno(i);
        cos = coss(i);
        N1 = bars(i,1);
        N2 = bars(i,2);
        A = bars(i,3);
        E = bars(i,4);
        I = bars(i,5);
        ro = bars(i,6);

        if strcmp(element_type,'plane frame')

            K11 = (E*A)/Ls;
            K22 = (12*E*I)/(Ls^3);
            K32 = (6*E*I)/(Ls^2);
            K23 = K32;
            K33 = (4*E*I)/Ls;
            K36 = (2*E*I)/Ls;
            K1 = [ K11,    0,     0, -K11,    0,     0;
                     0,  K22,   K23,    0, -K22,   K32;  
                     0,  K32,   K33,    0, -K32,   K36; 
                  -K11,    0,     0,  K11,    0,     0;
                     0, -K22,  -K23,    0,  K22,  -K32;
                     0,  K32,   K36,    0, -K32,   K33];

            M1 = ((ro*A*Ls)/420)*[140,      0,       0,  70,      0,       0;
                                    0,    156,   22*Ls,   0,     54,  -13*Ls;  
                                    0,  22*Ls,  4*Ls^2,   0,  13*Ls, -3*Ls^2; 
                                   70,      0,       0, 140,      0,       0;
                                    0,     54,   13*Ls,   0,    156,  -22*Ls;
                                    0, -13*Ls, -3*Ls^2,   0, -22*Ls,  4*Ls^2]; 

            Rot = [ cos,  sen,    0,    0,   0,   0;   % Rotation matrix
                    -sen,  cos,    0,    0,   0,   0;
                       0,    0,    1,    0,   0,   0;
                       0,    0,    0,  cos, sen,   0;
                       0,    0,    0, -sen, cos,   0;
                       0,    0,    0,    0,   0,   1];

            Klr = (Rot'*K1)*Rot;
            Mrr = (Rot'*M1)*Rot;

            gl1 = (3*N1)-2;
            gl2 = (3*N1)-1;
            gl3 = 3*N1;
            gl4 = (3*N2)-2;
            gl5 = (3*N2) -1;
            gl6 = 3*N2;

            K(gl1:gl3, gl1:gl3) = K(gl1:gl3, gl1:gl3) + Klr(1:3, 1:3);
            K(gl4:gl6, gl1:gl3) = K(gl4:gl6, gl1:gl3) + Klr(4:6, 1:3);
            K(gl1:gl3, gl4:gl6) = K(gl1:gl3, gl4:gl6) + Klr(1:3, 4:6);
            K(gl4:gl6, gl4:gl6) = K(gl4:gl6, gl4:gl6) + Klr(4:6, 4:6);

            M(gl1:gl3, gl1:gl3) = M(gl1:gl3, gl1:gl3) + Mrr(1:3, 1:3);
            M(gl4:gl6, gl1:gl3) = M(gl4:gl6, gl1:gl3) + Mrr(4:6, 1:3);
            M(gl1:gl3, gl4:gl6) = M(gl1:gl3, gl4:gl6) + Mrr(1:3, 4:6);
            M(gl4:gl6, gl4:gl6) = M(gl4:gl6, gl4:gl6) + Mrr(4:6, 4:6);

        elseif strcmp(element_type,'timoshenko beam')
            
            nu = 0.3;          % Poisson’s coefficient
            k = 0.5;           % Timoshenko’s shear factor               
            G = E/2/(1+nu);    % Elastic shear modulus

            fi = (12/Ls^2)*(E*I/(k*G*A));
            C = E*I/((1+fi)*Ls^3);

            K1 = C*[(A/I)*(1+fi),    0,           0, -(A/I)*(1+fi),     0,           0;
                               0,   12,        6*Ls,             0,   -12,        6*Ls;  
                               0, 6*Ls, (4+fi)*Ls^2,             0, -6*Ls, (2-fi)*Ls^2; 
                   -(A/I)*(1+fi),    0,           0,  (A/I)*(1+fi),     0,           0;
                               0,  -12,       -6*Ls,             0,    12,       -6*Ls;
                               0, 6*Ls, (2-fi)*Ls^2,             0, -6*Ls, (4+fi)*Ls^2];
                           
            M1 = ((ro*A*Ls)/(210*(1+fi)^2))*[70*(1+fi)^2,                           0,                           0, 35*(1+fi)^2,                           0,                           0;
                                                       0,       (70*fi^2)+(147*fi)+78, ((35*fi^2)+(77*fi)+44)*Ls/4,           0,        (35*fi^2)+(63*fi)+27,-((35*fi^2)+(63*fi)+26)*Ls/4;  
                                                       0, ((35*fi^2)+(77*fi)+44)*Ls/4, ((7*fi^2)+(14*fi)+8)*Ls^2/4,           0, ((35*fi^2)+(63*fi)+26)*Ls/4,-((7*fi^2)+(14*fi)+6)*Ls^2/4; 
                                             35*(1+fi)^2,                           0,                           0, 70*(1+fi)^2,                           0,                           0;
                                                       0,        (35*fi^2)+(63*fi)+27, ((35*fi^2)+(63*fi)+26)*Ls/4,           0,       (70*fi^2)+(147*fi)+78,-((35*fi^2)+(77*fi)+44)*Ls/4;
                                                       0,-((35*fi^2)+(63*fi)+26)*Ls/4,-((7*fi^2)+(14*fi)+6)*Ls^2/4,           0,-((35*fi^2)+(77*fi)+44)*Ls/4,((7*fi^2)+(14*fi)+8)*Ls^2/4]; 
            M2 = ((ro*I)/(30*(1+fi)^2*Ls))*[0,               0,                         0, 0,              0,                         0;
                                            0,              36,           -((15*fi)-3)*Ls, 0,            -36,           -((15*fi)-3)*Ls;
                                            0, -((15*fi)-3)*Ls, ((10*fi^2)+(5*fi)+4)*Ls^2, 0, ((15*fi)-3)*Ls,  ((5*fi^2)-(5*fi)-1)*Ls^2;
                                            0,               0,                         0, 0,              0,                        0;
                                            0,             -36,            ((15*fi)-3)*Ls, 0,             36,            ((15*fi)-3)*Ls;
                                            0, -((15*fi)-3)*Ls,  ((5*fi^2)-(5*fi)-1)*Ls^2, 0, ((15*fi)-3)*Ls, ((10*fi^2)+(5*fi)+4)*Ls^2];

            Mt = M1 + M2;

            Rot = [ cos,  sen,    0,    0,   0,   0;
                    -sen,  cos,    0,    0,   0,   0;
                       0,    0,    1,    0,   0,   0;
                       0,    0,    0,  cos, sen,   0;
                       0,    0,    0, -sen, cos,   0;
                       0,    0,    0,    0,   0,   1];

            Klr = (Rot'*K1)*Rot;
            Mrr = (Rot'*Mt)*Rot;

            gl1 = (3*N1)-2;
            gl2 = (3*N1)-1;
            gl3 = 3*N1;
            gl4 = (3*N2)-2;
            gl5 = (3*N2) -1;
            gl6 = 3*N2;

            K(gl1:gl3, gl1:gl3) = K(gl1:gl3, gl1:gl3) + Klr(1:3, 1:3);
            K(gl4:gl6, gl1:gl3) = K(gl4:gl6, gl1:gl3) + Klr(4:6, 1:3);
            K(gl1:gl3, gl4:gl6) = K(gl1:gl3, gl4:gl6) + Klr(1:3, 4:6);
            K(gl4:gl6, gl4:gl6) = K(gl4:gl6, gl4:gl6) + Klr(4:6, 4:6);
            
            M(gl1:gl3, gl1:gl3) = M(gl1:gl3, gl1:gl3) + Mrr(1:3, 1:3);
            M(gl4:gl6, gl1:gl3) = M(gl4:gl6, gl1:gl3) + Mrr(4:6, 1:3);
            M(gl1:gl3, gl4:gl6) = M(gl1:gl3, gl4:gl6) + Mrr(1:3, 4:6);
            M(gl4:gl6, gl4:gl6) = M(gl4:gl6, gl4:gl6) + Mrr(4:6, 4:6);

        elseif strcmp(element_type,'plane truss')
            
            Ke = (E*A)/Ls;
            K1 = Ke*[ 1, 0, -1, 0;
                      0, 0,  0, 0;  
                     -1, 0,  1, 0; 
                      0, 0,  0, 0];
                  
            Me = (ro*A*Ls)/6;
            M1 = Me*[ 2, 0,  1, 0;
                      0, 2,  0, 1;  
                      1, 0,  2, 0; 
                      0, 1,  0, 2];

            Rot = [ cos, sen,   0,   0;
                    -sen, cos,   0,   0;  
                       0, 0,   cos, sen; 
                       0, 0,  -sen, cos];

            Krr = (Rot'*K1)*Rot;
            Mrr = (Rot'*M1)*Rot;

            gl1 = (2*N1)-1;
            gl2 = 2*N1;
            gl3 = (2*N2)-1;
            gl4 = 2*N2;

            K(gl1:gl2, gl1:gl2) = K(gl1:gl2, gl1:gl2) + Krr(1:2, 1:2);
            K(gl3:gl4, gl1:gl2) = K(gl3:gl4, gl1:gl2) + Krr(3:4, 1:2);
            K(gl1:gl2, gl3:gl4) = K(gl1:gl2, gl3:gl4) + Krr(1:2, 3:4);
            K(gl3:gl4, gl3:gl4) = K(gl3:gl4, gl3:gl4) + Krr(3:4, 3:4);
            
            M(gl1:gl2, gl1:gl2) = M(gl1:gl2, gl1:gl2) + Mrr(1:2, 1:2);
            M(gl3:gl4, gl1:gl2) = M(gl3:gl4, gl1:gl2) + Mrr(3:4, 1:2);
            M(gl1:gl2, gl3:gl4) = M(gl1:gl2, gl3:gl4) + Mrr(1:2, 3:4);
            M(gl3:gl4, gl3:gl4) = M(gl3:gl4, gl3:gl4) + Mrr(3:4, 3:4);
            
        end
    end

    % Boundary conditions
    restr = zeros(1,dof*size(nodes,1));
    for i = 1:size(nodes,1)
      
        if strcmp(element_type,'timoshenko beam') || strcmp(element_type,'plane frame')

            RX = nodes(i,3);
            RY = nodes(i,4);
            RZ = nodes(i,5);

            gl1 = int32((3*i)-2);
            gl2 = int32((3*i)-1);
            gl3 = int32(3*i);

            if RX ==1
               restr(1,gl1)= 1;
            end

            if RY ==1
               restr(1,gl2)= 1;
            end

            if RZ ==1
               restr(1,gl3)= 1;
            end
            
        elseif strcmp(element_type,'plane truss')
 
            RX = nodes(i,3);
            RY = nodes(i,4);

            gl1 = int32((2*i)-1);
            gl2 = int32(2*i);

            if RX ==1
               restr(1,gl1)= 1;
            end

            if RY ==1
               restr(1,gl2)= 1;
            end
            
        end
        ccnt = restr==0;
        Kr = K(ccnt,ccnt);
        Mr = M(ccnt,ccnt);

    end

end
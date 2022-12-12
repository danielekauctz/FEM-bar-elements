% -------------------------------------------------------------------------
% Static Analysis: Nodal Displacements and Forces
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com
% -------------------------------------------------------------------------
% Input data:
% nodes: text file with information about structure nodes
    ... column 1 and 2: node X and Y coordinate
    ... column 3, 4 and 5: node boundary conditions (horizontal displacement, 
%             vertical displacement, rotation (if = 1 movement is restricted)
    ... column 6, 7 and 8: acting forces in the node (horizontal force, 
%             vertical force, bending moment
%       obs., truss elements have no rotation or bending moment, so there
%       is no columns 5 and 8
% element_type: bar element type
               ... 'plane frame'
               ... 'timoshenko beam'
               ... 'plane truss'
% K: stiffness matrix
% Kr: stiffness matrix with boundary conditions
% -------------------------------------------------------------------------
% Output data:
% disp: nodal displacements vector;
% forces: nodal forces vector;
% -------------------------------------------------------------------------
function [disp,forces] = static_analysis(nodes,element_type,K,Kr)

    if strcmp(element_type,'timoshenko beam') || strcmp(element_type,'plane frame')
        dof = 3;
        maxdof = dof*size(nodes,1);
        F = zeros(maxdof,1);
        restr = zeros(1,dof*size(nodes,1));
        
        for i = 1:size(nodes,1)
            RX = nodes(i,3);
            RY = nodes(i,4);
            RZ = nodes(i,5);
            Fx = nodes(i,6);
            Fy = nodes(i,7);
            Mz = nodes(i,8);

            gl1 = int32(3*i-2);
            gl2 = int32(3*i-1);
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
            
            if Fx ~= 0
                F(gl1) = Fx;
            end
            
            if Fy ~= 0
                F(gl2) = Fy;
            end
            
            if Fx ~= 0
                F(gl3) = Mz;
            end
        end
        
    elseif strcmp(element_type,'plane truss')
        dof = 2;
        maxdof = dof*size(nodes,1);
        F = zeros(maxdof,1);
        restr = zeros(1,dof*size(nodes,1));
        
        for i = 1:size(nodes,1)
            RX = nodes(i,3);
            RY = nodes(i,4);
            Fx = nodes(i,5);
            Fy = nodes(i,6);

            gl1 = int32((2*i)-1);
            gl2 = int32(2*i);
            
            if RX ==1
               restr(1,gl1)= 1;
            end

            if RY ==1
               restr(1,gl2)= 1;
            end
            
            if Fx ~= 0
                F(gl1) = Fx;
            end
            
            if Fy ~= 0
                F(gl2) = Fy;
            end
        end
        
    end
    ccnt = restr==0;
    
    Fr = F(ccnt);
        
    disp_r = inv(Kr)*Fr;
    disp = zeros(maxdof,1);
    disp(ccnt) = disp_r;

    forces = K*disp;

end
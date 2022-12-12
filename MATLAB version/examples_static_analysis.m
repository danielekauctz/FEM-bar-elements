%% FEM Bar Elements Problems: Static Analysis
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

clear; close; clc;
%% Example 1: Plane frame  --------------------------------------------------
% Input data
nodes = importdata('Input_planeframe_nodes.txt');
nan_nodes = find(isnan(nodes)); 
nodes(nan_nodes) = 0;

bars = importdata('Input_planeframe_bars.txt');
nan_bars = find(isnan(bars));
nodes(nan_bars) = 0;

element_type = 'plane frame';

% Figure
figure
plot_structure(nodes,bars,element_type);
hold off

% FEM Matrices
[K_pf,~,Kr_pf,~] = FEM_matrices(nodes,bars,element_type);

% Solution: Nodal Displacements and Forces
[disp_pf,forces_pf] = static_analysis(nodes,element_type,K_pf,Kr_pf);

%% Example 2: Truss structure ---------------------------------------------
% Input data
nodes = importdata('Input_truss_nodes.txt');
nan_nodes = find(isnan(nodes)); 
nodes(nan_nodes) = 0;

bars = importdata('Input_truss_bars.txt');
nan_bars = find(isnan(bars));
nodes(nan_bars) = 0;

element_type = 'plane truss';

% Figure
figure
plot_structure(nodes,bars,element_type);
hold off

% FEM Matrices
[K_truss,~,Kr_truss,~] = FEM_matrices(nodes,bars,element_type);

% Solution: Nodal Displacements and Forces
[disp_truss,forces_truss] = static_analysis(nodes,element_type,K_truss,Kr_truss);

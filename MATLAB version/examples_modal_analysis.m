%% FEM Bar Elements Problems: Modal Parameters
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

clear; close; clc;
%% Example 1: Timoshenko beam ---------------------------------------------
% Input data
nodes = importdata('Input_beam_nodes.txt');
nan_nodes = find(isnan(nodes)); 
nodes(nan_nodes) = 0;

bars = importdata('Input_beam_bars.txt');
nan_bars = find(isnan(bars));
nodes(nan_bars) = 0;

element_type = 'timoshenko beam';

% Figure
figure
plot_structure(nodes,bars,element_type);
hold off

% FEM Matrices
[~,~,Kr_beam,Mr_beam] = FEM_matrices(nodes,bars,element_type);

% Vibrational modes
[fn_beam,phi_beam] = modal_analysis(Kr_beam,Mr_beam);

%% Example 2: Plane frame  ------------------------------------------------
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
[~,~,Kr_pf,Mr_pf] = FEM_matrices(nodes,bars,element_type);

% Vibrational modes
[fn_pf,phi_pf] = modal_analysis(Kr_pf,Mr_pf);

%% Example 3: Truss structure ---------------------------------------------
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
[~,~,Kr_truss,Mr_truss] = FEM_matrices(nodes,bars,element_type);

% Vibrational modes
[fn_truss,phi_truss] = modal_analysis(Kr_truss,Mr_truss);
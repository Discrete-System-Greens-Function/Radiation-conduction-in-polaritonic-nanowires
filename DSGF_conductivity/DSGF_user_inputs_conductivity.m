%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSGF: User Inputs
%
%
% DESCRIPTION: This script should be edited with all user inputs for DSGF 
%              simulations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace and close all figures
clear, clc, close all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************START OF USER INPUTS***************************%

%*******************************DESCRIPTION*******************************%

% Short description of system you are modeling (this will be used to name
% the saved files)

description = 'test_nanowire';% nanowire spheres_test nanowire_nonuniform_test

%********************SELECTION OF TYPE OF SIMULATION *********************%

% Choose between sample or user_defined 

spatial_discretization_type = 'user_defined'; 
    
if strcmp('user_defined',spatial_discretization_type)
    
    %********************DISCRETIZATION OF THE SYSTEM*********************%
    %
    % Define the discretization for the system (2 group of objects).
    % The user should modify the discretization and delta_V parameters
    % according to the name of the file with the desired user-defined discretization. 
    % These files are generated using matlab scripts.
        
    discretization = "nanowire_cylinder_cut_12_Lchar_3.3e-08_d_0_cylinder_coating_N_12_Lchar_3.3e-08_N_3600_discretization"; %"nanowire_v1_200nm_N_1352_discretization"; % 
    delta_V = "nanowire_cylinder_cut_12_Lchar_3.3e-08_d_0_cylinder_coating_N_12_Lchar_3.3e-08_N_3600_delta_V_vector"; %"nanowire_v1_200nm_N_1352_delta_V_vector"; % 
    N_cross_section = 12; % 52 532 447
    structure = "coating"; 
    
    %{
    discretization = "nanowire_v1_200nm_N_1352_discretization"; % 
    delta_V = "nanowire_v1_200nm_N_1352_delta_V_vector"; % 
    N_cross_section = 52;
    %Choose between bare or coating 
    structure = "bare";
    %}
end

%********************************MATERIAL*********************************%
% Options:
%     'SiO2'
%     'SiC'
%     'Si3N4'
%     'user_defined'
if strcmp('bare',structure)
    material = Material.SiC;
elseif strcmp('coating',structure)
    material = Material.SiC_Au;
end    

%********************DIELECTRIC FUNCTION OF BACKGROUND********************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%************************FREQUENCY DISCRETIZATION*************************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)

% Insert the name of the file with the spectral discretization.
spectral_discretization = "SiC_fequencies_test.csv" ;
%spectral_discretization = "frequencies_all.csv" ;

%***********************TEMPERATURE OF EACH OBJECT************************%

%T = [400, 300]; % [K]
T = [300, 300]; % [K]


%****************TEMPERATURE FOR CONDUCTANCE CALCULATIONS*****************%

% Temperature at which the spectral conductance will be calculated.

%T_cond = [200,250,300,350,400]; % [K] 300;%
T_cond = [50,100,150,200,250,300,350,400]; % [K] 300;%


%*****************************DESIRED OUTPUTS*****************************%

% Save all Workspace variables in .mat file?
output.save_workspace = true;

% Save figures?
output.save_fig = true;

% figure format
output.figure_format = FigureFormat.fig;

% Output the total and spectral conductance for each bulk object?
output.conductance = true;

% Output the power dissipated in every subvolume?
output.power_dissipated_subvol = true;

% Output the power dissipated in each bulk object?
output.power_dissipated_bulk = true;

% Output the heatmap into slices?
output.heatmap_sliced = false;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = true;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = true;

                           
%***************************END OF USER INPUTS****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[omega] = readmatrix(append("Library/spectral_discretizations/", spectral_discretization));

if strcmp('sample',spatial_discretization_type) 
    delta_V = '';
    wave_type = "total";
elseif strcmp('user_defined',spatial_discretization_type)
    L_char = '';
    d = '';
    wave_type = "total";
end    


%DSGF_main(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, discretization_type, L_char, origin, delta_V);
%DSGF_main(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d);
%DSGF_main_conductivity(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d);
%DSGF_main_conductivity(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d, N_cross_section);

DSGF_main_conductivity_new(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d, N_cross_section, structure);

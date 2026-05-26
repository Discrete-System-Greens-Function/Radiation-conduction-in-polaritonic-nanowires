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

if strcmp('sample',spatial_discretization_type) 
    
    %********************DISCRETIZATION OF EACH OBJECT********************%
    % 
    % Define the discretization for each bulk object. In the sample, each bulk object needs 
    % its own discretization. The discretizations can be taken from the
    % pre-made samples or defined by the user.  
    %
    % The number at the end of the chosen discretization represents the 
    % number of subvolumes in that discretization. 
    %
    % Pre-made sample discretization options:
    %     Discretization.sphere_*
    %     Discretization.cube_*
    %
    % Example with two sample discretizations chosen:
    %      discretization = {Discretization.sphere_8, Discretization.sphere_8};
    %

    discretization = {Discretization.sphere_8, Discretization.sphere_8};
    
    %**************************SCALE EACH OBJECT**************************%

    % Characteristic length for scaling the discretized lattice of each bulk
    % object.
    %
    % If a pre-made sample is chosen, the characteristic length is:
    %     sphere: radius
    %     dipole: radius
    %     cube: side length
    %
    % If a user-defined input is chosen, the characteristic length is the
    % scaling factor of the user-input cubic lattice.
    %

    L_char = [50.e-9, 50.e-9]; % [m] %[50.e-9, 50.e-9]
   
    %**********************DISTANCE BETWEEN OBJECTS***********************%
    
    % Distance between the objects
    d =10.e-9; %[m]


    
elseif strcmp('user_defined',spatial_discretization_type)
    
    %********************DISCRETIZATION OF THE SYSTEM*********************%
    %
    % Define the discretization for the system (2 group of objects).
    % The user should modify the discretization and delta_V parameters
    % according to the name of the file with the desired user-defined discretization. 
    % These files are generated using matlab scripts.
        
    discretization = "nanowire_v1_200nm_N_1352_discretization"; % 
    delta_V = "nanowire_v1_200nm_N_1352_delta_V_vector"; % 
    N_cross_section = 52; %532 447
end

%********************************MATERIAL*********************************%
% Options:
%     'SiO2'
%     'SiC'
%     'Si3N4'
%     'user_defined'

material = Material.SiO2;
%material = Material.SiC;


%********************DIELECTRIC FUNCTION OF BACKGROUND********************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%************************FREQUENCY DISCRETIZATION*************************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)

% Insert the name of the file with the spectral discretization.
%spectral_discretization = "SiC_173_nonuniform_10_250_Trad_s.csv";
%spectral_discretization = "SiC_173_nonuniform_10_250_Trad_s.csv"; %
%spectral_discretization = "SiO2_100_uniform_5_25_um.csv" ;% 
%spectral_discretization = "SiO2_25_uniform_198_220_Trad_s.csv" ;%
%spectral_discretization = "SiO2_4_uniform_low.csv" ;
%spectral_discretization = "SiO2_120_nonuniform_5_25_um.csv" ;
%spectral_discretization = "SiO2_fequencies_f.csv" ;
spectral_discretization = "frequencies_all.csv" ;

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
output.DSGF_matrix = false;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = false;

                           
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
DSGF_main_conductivity(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d, N_cross_section);


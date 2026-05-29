function [ spectral_therm_conductivity, spectral_conductance ] = spectral_therm_conductivity_conductance_function( Trans, T_therm_cond, r, Lsub, omega, Nmin_SiC_1, Nmin_SiC_2, Nmax_SiC_1, Nmax_SiC_2, N_cross_section )
%function [ spectral_therm_conductivity ] = spectral_therm_conductivity_function( Trans, T_therm_cond, r, Lsub, omega )
% This function calculates the spectral thermal conductivity between two
% halves of a block of particles along each Cartesian direction.

% INPUTS:  Trans            (N x N) transmission coefficient matrix
%          T_therm_cond     Temperatuer at which to calculate thermal conductivity [K]
%          r                (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          d                Center-of-mass particle separation distance
%          omega            Radial frequency [rad/s]
%          ind_bulk_vector  Vector containing the index of the start of every bulk object in the r array of subvolume coordinates
%
% OUTPUTS: spectral_therm_conductivity   Spectral thermal conductivity [W/(K*m*(rad/s))]


% Constants
h_bar = 1.0545718e-34;      % Planck's constant [J*s]
k_b = 1.38064852e-23;       % Boltzmann constant [J/K]

% Total number of subvolumes in all bulk objects
[N,~] = size(r);

% % Cross sectional area 
%THE CROSS SECTIONAL AREA IS CALCULATED BASED ON THE NUMBER OF SUBVOLUMES
%IN THE CROSS-SECTIONAL AREA 
 A_x = N_cross_section*(Lsub)^2; 
% A_y = A_x;
% A_z = A_x;

% Summation along x-direction
summation_x = 0; % Initiate sum
summation_conductance = 0; % Initiate sum
for j = Nmin_SiC_1:Nmax_SiC_1 % Loop through subvolumes in volume B
    for i = Nmin_SiC_2:Nmax_SiC_2 % Loop through subvolumes in volume A        
        summation_x = summation_x + Trans(i,j)*(r(j,1) - r(i,1)); 
        summation_conductance = summation_x + Trans(i,j); %implemented on 05/29/2026
    end
end
%{
%compare
summation_xx = 0; % Initiate sum
for j = Nmin_SiC_2:Nmax_SiC_2 % Loop through subvolumes in volume B
    for i = Nmin_SiC_1:Nmax_SiC_1 % Loop through subvolumes in volume A        
        summation_xx = summation_xx + Trans(i,j)*(r(j,1) - r(i,1));        
    end
end
%}
%{

% Index particle locations
x_inds_A = find(r(:,1) > 0).';
x_inds_B = find(r(:,1) < 0).';
y_inds_A = find(r(:,2) > 0).';
y_inds_B = find(r(:,2) < 0).';
z_inds_A = find(r(:,3) > 0).';
z_inds_B = find(r(:,3) < 0).';
%{
% Number of particles along each direction
N_x = length(unique(r(:,1)));
N_y = length(unique(r(:,2)));
N_z = length(unique(r(:,3)));

% Cross sectional areas along each direction
A_x = N_y*N_z*(d^2);
A_y = N_x*N_z*(d^2);
A_z = N_x*N_y*(d^2);
%}

% % Cross sectional areas for a chain of spheres (Eric's solution)
% A_x = pi*(25e-9)^2;
% A_y = A_x;
% A_z = A_x;

x_inds_A = find(r(:,1) > 0).';
x_inds_B = find(r(:,1) < 0).';
% Summation along x-direction
summation_x = 0; % Initiate sum
for j = x_inds_B % Loop through subvolumes in volume B
    for i = x_inds_A % Loop through subvolumes in volume A        
        summation_x = summation_x + Trans(i,j)*(r(j,1) - r(i,1));        
    end
end
% For the current study, we are only calculating thermal conductivity along
x-direction.
% Summation along y-direction
summation_y = 0; % Initiate sum
for j = y_inds_B % Loop through subvolumes in volume B
    for i = y_inds_A % Loop through subvolumes in volume A        
        summation_y = summation_y + Trans(i,j)*(r(j,2) - r(i,2));        
    end
end

% Summation along z-direction
summation_z = 0; % Initiate sum
for j = z_inds_B % Loop through subvolumes in volume B
    for i = z_inds_A % Loop through subvolumes in volume A        
        summation_z = summation_z + Trans(i,j)*(r(j,3) - r(i,3));        
    end
end
%}

% Derivative of mean energy of an electromagnetic state with respect to
% temperature
num = ((h_bar*omega)^2).*exp((h_bar*omega)./(k_b.*T_therm_cond));
denom = k_b.*(T_therm_cond.^2).*((exp((h_bar*omega)./(k_b.*T_therm_cond)) - 1).^2);
dTheta_dT = num./denom;

spectral_therm_conductivity_x = -1.*(1/(A_x))*dTheta_dT*summation_x;%Livia removed 1/(2*pi)
% Consolidate into matrix
%spectral_therm_conductivity = [spectral_therm_conductivity_x, spectral_therm_conductivity_y, spectral_therm_conductivity_z];
spectral_therm_conductivity = [spectral_therm_conductivity_x]; %for the nanowire, we are only interested in the x-component

%%%%%%%%%%%%%% Version divided by 2 pi %%%%%%%%%%%%%%%%%
% Calculate spectral thermal conductivity [W/(K*m*(rad/s))]
%spectral_therm_conductivity_x = -1.*(1/(2*pi*A_x))*dTheta_dT*summation_x; % Lindsay version
%spectral_therm_conductivity_y = -1.*(1/(2*pi*A_y))*dTheta_dT*summation_y;
%spectral_therm_conductivity_z = -1.*(1/(2*pi*A_z))*dTheta_dT*summation_z;
%spectral_therm_conductivity_divided_by_2pi = [spectral_therm_conductivity_x]; %for the nanowire, we are only interested in the x-component
%%%%%%%%%%%%%% end Version divided by 2 pi %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Conductance %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spectral radiative conductance 
%To calculate conductance, it needs to be a bulk transmission coefficient 
%G_omega_bulk = dtheta_dT.*Trans_omega_bulk;
spectral_conductance = dTheta_dT.*summation_conductance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% end Conductance %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


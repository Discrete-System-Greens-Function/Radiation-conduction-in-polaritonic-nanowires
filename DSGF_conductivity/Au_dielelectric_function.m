function [ epsilon ] = Au_dielectric_function( omega )

% Dielectric function for free-electron gas (ideal) metal

%
%{
% example of range of frequencies
material = 'Au' ;
%ACS photonics SiC
min = 1.e14;
max = 5e14;
omega = linspace(min,max,10000);
%}

%extracted from DOI: 10.1007/s11468-015-0128-7

epsilon_0 = 9.84; 
omega_p = 13.71e15; % gold 
%tau = 1.e-14;
gamma = 1.09e14;%1/tau;

% Drude model of optical response of metals
epsilon = epsilon_0 - omega_p^2./(omega.^2+1i*omega.*gamma);


%{
% other set of parameters, I do not recall the source. 
omega_p = 13.71e15; % gold 
%tau = 1.e-14;
gamma = 4.05e13;%1/tau;

% Drude model of optical response of metals
epsilon = 1. - omega_p^2./(omega.^2+1i*omega.*gamma); 

%}

% case for large frequencies close to omega_p (not applicable to noble metals)
%epsilon = 1. - omega_p^2./(omega.^2); 

%epsilon_Au_real = real(epsilon_Au);
%epsilon_Au_imag = imag(epsilon_Au);
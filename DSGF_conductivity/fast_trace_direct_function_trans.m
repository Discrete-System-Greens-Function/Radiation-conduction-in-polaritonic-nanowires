function [ Trans] = direct_function_trans( omega, r, epsilon, epsilon_ref, delta_V, T_vector, ind_bulk, wave_type, T_conductance)

% This function calculates the system Green's function, the monochromatic 
% thermal power dissipated in each subvolume, and the monochromatic thermal
% power dissipated in each bulk object.

% INPUTS:  omega            Radial frequency [rad/s]
%          r                (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          epsilon          (N x 1) vector of all subvolume dielectric functions
%          epsilon_ref      Dielectric function for background reference medium (constant)
%          delta_V          (N x 1) vector of all subvolume sizes [m^3]
%          T_vector         (N x 1) vector of all subvolume temperatures [K]
%          ind_bulk         Indices of first subvolume in a given bulk object
%
% OUTPUTS: G_sys_2D         System Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)
%          Trans            Transmission coefficient matrix (3N x 3N)
%          Q_omega_bulk     Spectral heat dissipation within bulk objects [W/m]
%          Q_omega_subvol   Spectral heat dissipation within subvolumes [W/m*(m^3)]
%

% Determine total number of subvolumes
[N,~] = size(r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate background medium Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[ ~, G_0_2D] = G_0_function(r, omega, delta_V, epsilon_ref, N, wave_type); % original
[~, G_0_2D] = fast_trace_G_0_function(r, omega, delta_V, epsilon_ref, N, wave_type); % modification: parfor
                                                                    % modification for triangular matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate deterministic interaction matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A_2D, ~ ] = A_matrix_function(G_0_2D, omega, epsilon, epsilon_ref, delta_V, N);

%[ A_packed, ~ ] = fast_trace_A_matrix_function(G_0_4D, omega, epsilon, epsilon_ref, delta_V, N);
%display_memory_consumption(struct2cell(whos));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function (Original) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = toc;
%G_sys_2D = A_2D\G_0_2D;  % System Green's function matrix

G_sys_2D = A_2D\G_0_2D(:, 1:(3*N/2));  % System Green's function matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function using GMRES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computes half of SGF (the ones that are necessary of the trace) saving some memory. 
%The tolerance needed to be increased, unclear if we need to keep increasing the tolerance as we increase N. 
%G_0_3N_1 = ones(3*N,1);
%i=0;
%x0 = zeros(size(G_0_2D,1),1); % first guess


%{
%use of the pre-conditioner increased the memory usage
setup.type = 'ilutp'; % ILU with threshold pivoting
A = sparse(A_2D);
clear A_2D;
[L,U] = ilu(A, setup);
%}
%M = diag(diag(A_2D));  % diagonal preconditioner


%{
for column = 1:3*N/2    
    %G_0_3N_1(:,1) = G_0_2D(:,column);
    %G_0_3N_1(:,1) = G_0_2D((3*N/2)+1:3*N,column);

    % Vectorize B
    %b_vec = G_0_3N_1(:);
    %b_vec = G_0_2D(:,column);
      
    % Solve using GMRES
    tol = 1e-8; % 1e-6  1e-10
    maxit = 3*N; % max iterations
    %[x_vec, flag] = gmres(A_2D, b_vec, [], tol, maxit); %
    [x_vec, flag] = gmres(A_2D, G_0_2D(:,column), [], tol, maxit); %
    
    %[x_vec,flag] = gmres(A_2D, G_0_2D(:,column), [], tol, maxit, M); %takes more time and more memory
    
    %[x_vec, flag] = gmres(A, G_0_2D(:,column), [], tol, maxit, L, U); %increased memory usage
    
    %[x_vec, flag] = gmres(A, G_0_2D(:,column), [], tol, maxit,L,U, x0); %
    %x0 = x_vec; % warm start for next RHS
    
    %display_memory_consumption(struct2cell(whos));
% Reshape solution back to N x N
%X = reshape(x_vec, N, N);
%G_sys_2D = reshape(x_vec, 3*N, column);
    G_sys_2D(:,column) = x_vec;
    %display_memory_consumption(struct2cell(whos));
    %i = i+1;
end
%display(i);
%}


%display_memory_consumption(struct2cell(whos));

%clear A;
%clear L;
%clear U;
%clear M;


clear G_0_2D;
clear A_2D;


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function using GMRES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_0_3N_1 = ones(3*N,1);
for column = 1:3*N    
    G_0_3N_1(:,1) = G_0_2D(:,column);

    % Vectorize B
    b_vec = G_0_3N_1(:);
    % Solve using GMRES
    tol = 1e-6; %1e-10
    maxit = 3*N; % max iterations
    [x_vec, flag] = gmres(A_2D, b_vec, [], tol, maxit); %
    %display_memory_consumption(struct2cell(whos));
% Reshape solution back to N x N
%X = reshape(x_vec, N, N);
%G_sys_2D = reshape(x_vec, 3*N, column);
    G_sys_2D(:,column) = x_vec;
    %display_memory_consumption(struct2cell(whos));

end

clear G_0_3N_1;
clear G_0_2D;
clear A_2D;
%display_memory_consumption(struct2cell(whos));
%}

t2 = toc;

disp(['Time for matrix inversion = ' num2str(t2-t1) ' s = ' num2str((t2-t1)/60) ' minutes'])

% Determine total number of subvolumes
[N,~] = size(r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function using GMRES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
t1 = toc;

N_reshape = size(A_2D,1);
% Vectorize B
b_vec = G_0_2D(:);
clear G_0_2D;

% Function that applies K*x without building K
Af = @(x) reshape(A_2D * reshape(x, N_reshape, N_reshape), N_reshape^2, 1);

% --- Preconditioner ---
if issparse(A_2D)
    % For sparse A, use incomplete LU
    setup.type = 'ilutp';
    setup.droptol = 1e-3;   % tune droptol for your problem
    [L,U] = ilu(sparse(A_2D), setup);
    Msolve = @(x) reshape(U \ (L \ reshape(x, 3*N, 3*N)), 9*N^2, 1);
else
    % For dense A, use full LU
    [L,U,P] = lu(A_2D);
    Msolve = @(x) reshape(U \ (L \ (P * reshape(x, 3*N, 3*N))), 9*N^2, 1);
end
clear A_2D;
%display_memory_consumption(struct2cell(whos));

% --- GMRES solve ---
tol = 1e-8;
maxit = 9*N^2; %200;        % max total iterations
restart = 50;       % restart size
[x_vec, flag, relres, iter] = gmres(Af, b_vec, restart, tol, maxit, Msolve);

%display_memory_consumption(struct2cell(whos));

clear Af;
clear b_vec;

% Reshape solution back to N x N
G_sys_2D = reshape(x_vec, 3*N, 3*N); % Model: X = reshape(x_vec, N, N);
clear x_vec;

t2 = toc;

disp(['Time for matrix inversion = ' num2str(t2-t1) ' s = ' num2str((t2-t1)/60) ' minutes'])

% Determine total number of subvolumes
[N,~] = size(r);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transmission coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t3 = toc;
%[ Trans ] = trans_coeff_function(G_sys_2D, omega, epsilon, delta_V, N);

%key difference: calculates transmission coefficient of half the terms
[ Trans ] = fast_trace_trans_coeff_function(G_sys_2D, omega, epsilon, delta_V, N);
%display_memory_consumption(struct2cell(whos));
t4 = toc;



disp(['Time to calculate transmission coefficient = ' num2str(t4-t3) ' s = ' num2str((t4-t3)/60) ' minutes'])







function [ G_sys_2D, Trans, Q_omega_bulk, Q_omega_subvol ] = direct_function( omega, r, epsilon, epsilon_ref, delta_V, T_vector, ind_bulk, wave_type )

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

[ ~, G_0_2D] = G_0_function(r, omega, delta_V, epsilon_ref, N, wave_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate deterministic interaction matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A_2D, ~ ] = A_matrix_function(G_0_2D, omega, epsilon, epsilon_ref, delta_V, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = toc;
G_sys_2D = A_2D\G_0_2D;  % System Green's function matrix
display_memory_consumption(struct2cell(whos));
%{
G_0_3N_1 = ones(3*N,1);
for column = 1:3*N    
    G_0_3N_1(:,1) = G_0_2D(:,column);

    % Vectorize B
    b_vec = G_0_3N_1(:);
    % Solve using GMRES
    tol = 1e-6; %1e-10
    maxit = 3*N; % max iterations
    [x_vec, flag] = gmres(A_2D, b_vec, [], tol, maxit); %
    display_memory_consumption(struct2cell(whos));
% Reshape solution back to N x N
%X = reshape(x_vec, N, N);
%G_sys_2D = reshape(x_vec, 3*N, column);
    G_sys_2D(:,column) = x_vec;
    %display_memory_consumption(struct2cell(whos));

end
%}
clear G_0_2D;
clear A_2D;

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


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First attempt using GMRES                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using GMRES without preconditioners
% Build the Kronecker product system
I_N = speye(3*N);  % identity matrix
K = kron(I_N, A_2D); % (I ⊗ A)

% Vectorize B
b_vec = G_0_2D(:);


% Solve using GMRES
tol = 1e-10;
maxit = 9*N^2; % max iterations
[x_vec, flag] = gmres(K, b_vec, [], tol, maxit); %

% Reshape solution back to N x N
%X = reshape(x_vec, N, N);
G_sys_2D = reshape(x_vec, 3*N, 3*N);

% Check residual
%residual = norm(A*X - B);
residual = norm(A_2D*G_sys_2D - G_0_2D);
disp(['Residual: ', num2str(residual)]);
%}

t2 = toc;

disp(['Time for matrix inversion = ' num2str(t2-t1) ' s = ' num2str((t2-t1)/60) ' minutes'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transmission coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t3 = toc;
[ Trans ] = trans_coeff_function(G_sys_2D, omega, epsilon, delta_V, N);
display_memory_consumption(struct2cell(whos));
t4 = toc;

disp(['Time to calculate transmission coefficient = ' num2str(t4-t3) ' s = ' num2str((t4-t3)/60) ' minutes'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate monochromatic power dissipated in bulk objects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Q_omega_bulk, Q_omega_subvol ] = Q_function(Trans, omega, T_vector, N, ind_bulk);




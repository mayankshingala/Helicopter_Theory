close all
clear all
clc

data;
pilot_inputs;
gaussian;

%% function to solve for β**, β*, β (Newmark's OR Harmonic balance)

% function to solve Harmonic balance
[beta_0_HB,beta_1c_HB, beta_1s_HB] = harmonic_balance(gamma,nu_beta,mu,theta_0,theta_tw,theta_1s,theta_1c,lambda);

beta_0 = beta_0_HB;
beta_1c = beta_1c_HB;
beta_1s = beta_1s_HB;

% To calculate β**, β*, β
n = 1;
for psi = (0:1:360)*(pi/180)        % for step of 1 deg
    beta(n) = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);
    beta_star(n) = -beta_1c*sin(psi) + beta_1s*cos(psi);
    beta_dstar(n) = -beta_1c*cos(psi) - beta_1s*sin(psi);
    beta_dot(n) = omega * beta_star(n);
    beta_ddot(n) = (omega^2) * beta_dstar(n);
    n = n + 1;
end

%% Fixed Point Iteration to calculate uniform inflow (λ)
CT_hover  = Weight / (rho*A*((omega*R)^2));             % Calculation of Coefficient of thrust (at hover)
CT = CT_hover;
% lambda_0 = sqrt(CT/2);
lambda = inflow(mu,alpha_s,CT);
% lambda = inflow(lambda_0,mu,alpha_s,CT);

%% Assigning Ψ values for all 4 blades over one revolution
psi_r = 0 : 1*(pi/180) : 2*pi;              % Vector of all azimuthal locations (Ψ) of blades with step of 1 deg

% Matrix of all azimuthal locations (Ψ) of blades...
% (Here 4 bladed rotor are separated by 90°)
psi_all = [psi_r;
           psi_r + (pi/2); 
           psi_r + (pi); 
           psi_r + (3*pi/2)];       % Create general case for Nb bladed rotor

%% Empty matrices to collect vaues of fx,fy,fz,mx,my,mz for rotating & fixed frame
[row, col] = size(psi_all);      % row- No. of rows, col-No. of columns
fx = zeros(row,col-1);
fy = zeros(row,col-1);
fz = zeros(row,col-1);

mx = zeros(row,col-1);
my = zeros(row,col-1);
mz = zeros(row,col-1);

Fx_m_fixed = zeros(row,col-1);
Fy_m_fixed = zeros(row,col-1);
Fz_m_fixed = zeros(row,col-1);

Mx_m_fixed = zeros(row,col-1);
My_m_fixed = zeros(row,col-1);
Mz_m_fixed = zeros(row,col-1);

% Calculation of m (mass/length)
I_beta = (rho*Cla*C*R^4)/gamma;
m = (3*I_beta) / ((R-e)^3);

integration;

plotting;
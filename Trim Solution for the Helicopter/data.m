%% Given data

D = 1.225;                  % Density of air [kg/m^3]
Nb = 4;                     % Number of blades
R = 8.18;                   % Blade radious [m]
C = 0.46;                   % Blade chord [m]
Cd0 = 0.01;                 % Profile drag coefficient
Cla = 5.73;                 % Lift curve slope
Omega = 27*pi;              % Rotor angular speed [rad/s]
nu_beta = 1.04;             % Blade flap frequency [/rev]
% nu_beta = nu_beta/(2*pi);   % Blade flap frequency [/rad]
gamma = 8;                  % Lock number
Weight = 70000;             % Weight of ACFT [N]
theta_tw = 0;               % blade twist rate [/rad OR /deg]

R_cut = 0;                  % Root cut-out (Considered for general case)
%% Given data

rho = 1.225;                % Density of air [kg/m^3]
Nb = 4;                     % Number of blades
R = 8.18;                   % Blade radious [m]
C = 0.46;                   % Blade chord [m]
Cd0 = 0.01;                 % Profile drag coefficient
Cla = 5.73;                 % Lift curve slope
omega = 27;                 % Rotor angular speed [rad/s]
nu_beta = 1.04;             % Blade flap frequency [/rev]
gamma = 8;                  % Lock number
Weight = 70000;             % Weight of ACFT [N]
theta_tw = 0;               % blade twist rate [/rad OR /deg]

% for question-02
% h = 1.83;                   % Vertical separation between Hub plane and CG [m]
% lT = 9.75;                  % Distance between CG and tail rotorhub center [m]
% f = 1.85;
% x_cg = -0.6;                % CG location w.r.t. hub in Longitudinal direction
% y_cg = 0;                   % CG location w.r.t. hub in Lateral direction
% Mx_F = 0;                   % Pitching moment due to fuselage w.r.t. longitudinal axis (x)
% My_F = 0;                   % Pitching moment due to fuselage w.r.t. lateral axis (y)
% k_h = 1.15;                 % Induced power factor (for hover)

%
R_cut = 0.06*R;             % Root cut-out (Considered 6% from data)
e = R_cut/R;                % Non-dimentional root cut-out
g = 9.81;                   % Gravitational acceleration [m/s^2]

Cd = Cd0;                   % Assumption
A = pi*R^2;                 % Rotor disk area
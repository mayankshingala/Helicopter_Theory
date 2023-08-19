%% Given data
% Constant values
d = 1.225;      % Density of air, ρ [kg/m3]
Nb = 4;         % Number of blades
R = 6;          % Blade Radius [m]
c = 0.4;        % Blade Chord [m]
Cd_0 = 0.01;    % Profile Drag coefficient
Cl_a = 5.73;    % Lift curve slope, Clα
omg = 10*pi;    % Rotor angular rate, Ω [rad/sec]
CT_h = 0.006;   % Hover thrust coefficient
R_cut = 0.2*R;  % Root cut-out
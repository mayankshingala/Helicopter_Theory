close all
clear all
clc

data;
inputs;
gaussian;
beta_p = 0;             % Pre-cone angle is 0

%% Calculated data
rev = 5;                    % Number of revolutions
d_psi = [deg2rad(2) deg2rad(3) deg2rad(5) deg2rad(10) deg2rad(15)];         % Time steps given in problem

% Newmark's Parameters
beta_nm = 1/4;
gamma_nm = 1/2;

% Harmonic Balance(HB) coefficients 
% (P*R = Q) Solved using back slash operator;
beta_0_HB = (gamma/(nu_beta^2)) * ( theta_0*(1+mu^2)/8 + (theta_tw/10)*(1 + (5/6)*(mu^2)) + mu*theta_1s/6 - lambda/6) ;

P = [1 (1+(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1))); 
    -(1-(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1))) 1];

Q = (gamma/(nu_beta^2 - 1))*[ (theta_1c/8)*(1+(mu^2)/2) - (mu/6)*beta_0_HB;
    (theta_1s/8)*(1-(mu^2)/2) + (mu/3)*theta_0 - (mu/4)*lambda + (mu^2/4)*theta_1s + mu*theta_tw/4];

R = P\Q;

beta_1c_HB = R(1);
beta_1s_HB = R(2);

%% Variable arrays declarations
M_beta = zeros(rev*2*pi/d_psi(1) + 1, 5);
beta = zeros(rev*2*pi/d_psi(1) + 1, 5);           % includes initial condition for Newmark's
beta_1 = zeros(rev*2*pi/d_psi(1) + 1, 5);
beta_2 = zeros(rev*2*pi/d_psi(1) + 1, 5);
beta_HB = zeros(rev*2*pi/d_psi(1) + 1, 5);
% beta(1,:) = 0.1774;                               % harmonic balance initial condition for better solution
beta(1,:) = 0;
%% Numerical Integration
k = 0;
for j = 1:1:5
    for psi = 0 : d_psi(j) : (rev*2*pi-d_psi(j))      % defined over whole azimuth (ψ) by number of revolutions
        k = k+1;

        % Newmark's Algorithm
        % M_BETA_bar (Mβ)
        M_beta(k,j) = ( (1/8) + (mu/3)*sin(psi) + ((mu^2)/4)*(sin(psi)^2))*...
                      (theta_0 +...
            theta_1c*cos(psi)+theta_1s*sin(psi))+...
            theta_tw*((1/10)+(mu^2)/6*(sin(psi)^2)+mu/4*sin(psi))...
            - lambda*((1/6)+mu/4*sin(psi)) - ...
            beta_1(k,j)*((1/8)+(mu/6)*sin(psi))...
            - mu*beta(k,j)*cos(psi)*((1/6)+(mu/4)*sin(psi));

        % β**
        if k == 1
            beta_2(k,j) = gamma*M_beta(k,j) - (nu_beta^2)*beta(k,j);
        end

        % Expression for β**(n+1), β*(n+1), β(n+1) which uses β**(n) and β*(n)
        beta_2(k+1,j) = (gamma*M_beta(k,j) -...
                        (nu_beta^2)*(d_psi(j)*beta_1(k,j) + beta(k,j) + ((d_psi(j)^2)/2)*(1 - 2*beta_nm)*beta_2(k,j)))/...
                        ( 1 + (nu_beta^2)*(d_psi(j)^2)*beta_nm );

        beta_1(k+1,j) = beta_1(k,j) + ...
                        d_psi(j)*( (1 - gamma_nm)*beta_2(k,j) + gamma_nm*beta_2(k+1,j) );

        beta(k+1,j) = beta(k,j) + ...
                      d_psi(j)*beta_1(k,j) +...
                      ((d_psi(j)^2)/2)*( (1 - 2*beta_nm)*beta_2(k,j) + 2*beta_nm*beta_2(k+1,j) );

        % Harmonic Balance
        beta_HB(k,j) = beta_0_HB + beta_1c_HB*cos(psi) + beta_1s_HB*sin(psi);
    end
    k = 0;
    
    % Plots
%     hold on
    psi = [0 : d_psi(j) : rev*2*pi-d_psi(j)];
    plot(psi, beta(1:length(psi),j))
    hold on
end

plot_2;

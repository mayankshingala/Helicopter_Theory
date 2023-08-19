close all
clear all
clc

data;
inputs;
gaussian;

%% calculated data
n = 20;                             % Radial stations
r_cut = (R_cut/R);                  % Non-dimensionalized Root cut-out
ds = (1-r_cut)/n;                   % Non-dimentional length of each segments
rev = 5;                            % Number of revolutions graph is defined for
d_psi = 5*pi/180;                   % Time steps
m = 1;                              % Azimuth angle initializer

%% Variables declarations to collect calculated data
f_M_beta = zeros(1,6);
M_beta = zeros(1, rev*2*pi/d_psi + 1);
M_beta_exact = zeros(1, rev*2*pi/d_psi + 1);
beta = zeros(1, rev*2*pi/d_psi + 1);
beta_star = zeros(1, rev*2*pi/d_psi + 1);

%% Numerical Integration
for psi = 0 : d_psi : rev*2*pi              % defined over whole azimuth (ψ)
    for j = 1:1:n

        a = r_cut + ds*(j-1);               % Lower limit of integration
        b = a + ds;                         % Upper limit of integration
        t = ((b-a)/2)*Node + ((b+a)/2);     % change of variables: new locations are t(i) from x(i)

        theta = theta_0 + theta_tw*(t - 0.75) + theta_1s*sin(psi) + theta_1c*cos(psi);      % Pilot input
        beta = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);                                % Response of blade (in pitch)
        
        beta_star = -beta_1c*sin(psi) + beta_1s*cos(psi);       % β* for exact solution

        for i = 1:1:6               % 6 point Guassian quadrature
            f_M_beta(i) = (1/2)*t(i) * ((t(i) + mu.*sin(psi)).^2 * theta(i) ...
                          - (lambda*cos(beta) + t(i)*beta_star + mu*sin(beta).*cos(psi))*(t(i) + mu.*sin(psi)));
            M_beta(m) = M_beta(m) + W(i)*f_M_beta(i)*((b-a)/2);
        end

    end

%     M_beta_exact(m) = M_beta (mu,psi,theta_0,theta_1c,theta_1s,theta_tw,lambda,beta_star)
    M_beta_exact(m) = ( (1/8) + (mu/3)*sin(psi) + (((mu^2)*(sin(psi))^2)/4) )* ...
                      ( theta_0 + theta_1c*cos(psi) + theta_1s*sin(psi) ) + ...
                      theta_tw * ( (1/10) + ((mu^2)/6)*((sin(psi))^2) + ((mu/4)*sin(psi)) ) - ...
                      lambda*( (1/6) + (mu/4)*sin(psi) ) - ...
                      beta_star*( (1/8) + (mu/6)*sin(psi) ) - ...
                      mu*beta*cos(psi)*( (1/6) + (mu/4)*sin(psi) );
    m = m + 1;
end

plot_1;

function [beta_0_HB,beta_1c_HB, beta_1s_HB] = harmonic_balance(gamma,nu_beta,mu,theta_0,theta_tw,theta_1s,theta_1c,lambda)

% Harmonic Balance(HB) coefficients 
% (P*R = Q) Solved using back slash operator;

beta_0_HB = (gamma/(nu_beta^2)) * ( theta_0*(1+mu^2)/8 + (theta_tw/10)*(1 + (5/6)*(mu^2)) + mu*theta_1s/6 - lambda/6) ;

A11 = 1;
A12 = (1+(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1)));
A21 = -(1-(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1)));
A22 = 1;

P = [A11, A12; A21, A22];

B11 = (theta_1c/8)*(1+(mu^2)/2) - (mu/6)*beta_0_HB;
B12 = (theta_1s/8)*(1-(mu^2)/2) + (mu/3)*theta_0 - (mu/4)*lambda + (mu^2/4)*theta_1s + mu*theta_tw/4;

Q = (gamma/(nu_beta^2 - 1)) * [B11; B12];

% P = [1 (1+(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1))); 
%     -(1-(mu^2)/2)*(gamma/(8*(nu_beta^2 - 1))) 1];

% Q = (gamma/(nu_beta^2 - 1))*[ (theta_1c/8)*(1+(mu^2)/2) - (mu/6)*beta_0_HB;
%     (theta_1s/8)*(1-(mu^2)/2) + (mu/3)*theta_0 - (mu/4)*lambda + (mu^2/4)*theta_1s + mu*theta_tw/4];

R = P\Q;

beta_1c_HB = R(1);
beta_1s_HB = R(2);

end
clear all
close all
clc

helidata_BEMT
n = 100;          % number of radial segments

%% calculated data

r_ct = 0.2;                                             % Non-dimensional root cutout
sgm = Nb*c*(R - R_cut)/(pi*R*R);                        % Rotor Solidity
lambda_h = sqrt(CT_h/2);                                % Induced inflow for hover
dy= (1-0.2)/n;                                          % Non-dimensional length of each segment
theta_0 = ((6*CT_h)/(sgm*Cl_a)) + (3/2)*sqrt(CT_h/2);   % Reference pitch angle from BET & MT
error = 10^-6;                                          % Error acceptance value

gaussian_BEMT

%% Numerical Solution

% Initialization of values
k = 1;                           % Array initializer for plots
CT(n) = 0;
CPi(n) = 0;
f_ct(6) = 0;
f_cp(6) = 0;

CT_TL(n) = 0;                    % While considering tip loss effect 
CPi_TL(n) = 0;
f_ct_TL(6) = 0;
f_cp_TL(6) = 0;

lambda_collect = 0;              % inflow analytical
lambda_collect_TL = 0;           % inflow using Prandtl's function

for i = 1:1:n           % Blade span is divided into n number of radial segments

    % Gaussian Quadrature for each segment
    a = 0.2 + dy*(i-1);                         % Integration lower limit
    b = 0.2 + dy*i;                             % Integration upper limit
    t = ((b-a)/2)*Node + ((b+a)/2);             % Change of variables, New locations are t(i)
    theta_tw = 0;                               % For no twist case
    theta = theta_0 + theta_tw*(t - 0.75);      % Final θ for initial guess

    for j = 1:1:6
        % No twist: θ_tw = 0 deg
        theta_tip = theta(j)*t(j);                                          % θ_tip calculation
        lambda_t = (sgm*Cl_a/16)*(sqrt(1+(32*theta_tip/(sgm*Cl_a))) - 1);   % inflow variation with r
        e = lambda_h;                                                       % error initialization
        lambda_TL_t = lambda_t;                                             % Initially starting with F=1 (No tip loss);

        while e > error
                f = (Nb/2)*(1-t(j))/lambda_TL_t;
                F = (2/pi)*acos(exp(-f));                                                   % Tip loss coefficients
                lambda_clc = (sgm*Cl_a/(16*F))*(sqrt(1+(32*F*theta_tip/(sgm*Cl_a))) - 1);   % Prandtl's tip loss function
                e = abs(lambda_TL_t - lambda_clc);                                          % Error calculation
                lambda_TL_t = lambda_clc;                                                   % λ for next iteration
        end

        lambda_collect = lambda_collect + lambda_t;                      % Collect λ for plotting
        lambda_collect_TL = lambda_collect_TL + lambda_TL_t;
        dCT_dr(j) = 0.5*sgm*Cl_a*(theta_tip - lambda_t)*t(j);            % CT distribution
        dCT_dr_TL(j) = 0.5*sgm*Cl_a*(theta_tip - lambda_TL_t)*t(j);
        
        f_ct(j) = dCT_dr(j);
        f_ct_TL(j) = dCT_dr_TL(j);
        
        CT(i) = CT(i) + W(j)*f_ct(j)*((b-a)/2);             % performing numerical integration for CT
        CT_TL(i) = CT_TL(i) + W(j)*f_ct_TL(j)*((b-a)/2);
        
        f_cp(j) = (dCT_dr(j)*lambda_t) + (0.5*sgm*Cd_0*(t(j)^3));
        f_cp_TL(j) = (dCT_dr_TL(j)*lambda_TL_t) + (0.5*sgm*Cd_0*(t(j)^3));
%         f_cp(j) = (0.5*sgm*Cl_a*(theta_tip - lambda_t)*t(j)*lambda_t) + (0.5*sgm*Cd_0*(t(j)^3));
%         f_cp_TL(j) = (0.5*sgm*Cl_a*(theta_tip - lambda_TL_t)*t(j)*lambda_TL_t) + (0.5*sgm*Cd_0*(t(j)^3));
        
        CPi(i) = CPi(i) + W(j)*f_cp(j)*0.5*(b-a); % performing numerical integration for CPi
        CPi_TL(i) = CPi_TL(i) + W(j)*f_cp_TL(j)*0.5*(b-a);
    end
    
    loc(k) = t(3);                            % Locations for plot
    lambda_mean(k) = lambda_collect/6;        % Inflow analytical  
    lambda_mean_TL(k) = lambda_collect_TL/6;  % Inflow using Prandtl's function
    
    lambda_collect = 0;                       % collect values to zero for next iteration
    lambda_collect_TL = 0;
    
    % The plot variables for CT, CPi for both cases
    CT_mean(k) = mean(dCT_dr);
    CT_TL_mean(k) = mean(dCT_dr_TL);
    
    CP_mean(k) = CPi(i);
    CP_TL_mean(k) = CPi_TL(i);
    k = k+1;                                  % Incrementing the location
end

results_BEMT_Q2
plot_BEMT_Q2


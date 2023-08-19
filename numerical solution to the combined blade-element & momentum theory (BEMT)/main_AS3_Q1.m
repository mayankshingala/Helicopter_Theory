clear all
close all
clc

helidata_BEMT
n = 10;         % number of radial segments

%% Calculated data

R_cut = 0.2;                                            % Non-dimensional root cutout
sgm = Nb*c*(R-R_cut)/(pi*R*R);                          % Rotor Solidity
lambda_h = sqrt(CT_h/2);                                % Induced inflow for hover
dy = (1-0.2)/n;                                         % Non-dimensional length of each segment
theta_0 = ((6*CT_h)/(sgm*Cl_a)) + (3/2)*sqrt(CT_h/2);   % Reference pitch angle from BET & MT
theta_75 = theta_0;                                     % Reference pitch angle considered at 75% of R
theta_tip = (4*CT_h/sgm/Cl_a) + sqrt(CT_h/2);           % Ideal Twist solution for pitch angle
c_tip = 0.5*c;                                          % Ideal Taper assumption
sgm_tip = Nb*c_tip*(R - R_cut)/(pi*R*R);                % Rotor solidity for ideal taper

gaussian_BEMT

%% Numerical Values Initialization

k = 1;                  % Array initializer for mean values

CT_1 = zeros(1,n);      % Initialisation for thrust coefficients
CT_2 = zeros(1,n);
CT_3 = zeros(1,n);
CT_4 = zeros(1,n);

CPi_1 = zeros(1,n);     % Initialisation for induced power coefficient (or torque coefficient)
CPi_2 = zeros(1,n);
CPi_3 = zeros(1,n);
CPi_4 = zeros(1,n);

CPp_1 = zeros(1,n);     % Initialisation for profile power coefficient (or torque coefficient)
CPp_2 = zeros(1,n);
CPp_3 = zeros(1,n);
CPp_4 = zeros(1,n);

f_ct_1 = zeros(1,6);    % Gaussian function initialisation for CT (Thrust coefficient)
f_ct_2 = zeros(1,6);
f_ct_3 = zeros(1,6);
f_ct_4 = zeros(1,6);

f_cp_1 = zeros(1,6);    % Gaussian function initialisation for CPi (Induced power coefficient) 
f_cp_2 = zeros(1,6);
f_cp_3 = zeros(1,6);
f_cp_4 = zeros(1,6);

f_cp_p1 = zeros(1,6);   % Gaussian function initialisation for CPp (Profile power coefficient)
f_cp_p2 = zeros(1,6);
f_cp_p3 = zeros(1,6);
f_cp_p4 = zeros(1,6);

lambda_collect_1 = 0;
lambda_collect_2 = 0;
lambda_collect_3 = 0;
lambda_collect_4 = 0;

for i = 1:1:n           % Blade span is divided into n number of radial segments

    % Gaussian Quadrature for each segment
    a = 0.2 + dy*(i-1);                             % Lower limit of integration
    b = 0.2 + dy*i;                                 % Upper limit of integration
    t = ((b-a)/2)*Node + ((b+a)/2);                 % Change of variables, New locations are t(i)
    theta_tw_2 = 0;                                 % Q1-(ii), θ_tw (twist rate) = 0deg
    theta_2 = theta_75 + theta_tw_2*(t - 0.75);     % θ for θ_tw = 0
    theta_tw_3 = -15*pi/180;                        % Q1-(ii), θ_tw (twist rate) = -15deg
    theta_3 = theta_75 + theta_tw_3*(t - 0.75);     % θ for θ_tw = -15deg
    sgm_4 = (Nb*(c_tip./t)*(R - R_cut))/(pi*R*R);   % Q1-(iv), Ideal twist & Ideal taper

    for j=1:1:6

        % Case-1: Ideal Twist: θ_tip

        theta_tip_1 = theta_tip;                                                % θ_tip calculation
        lambda_t_1 = (sgm*Cl_a/16)*(sqrt(1+(32*theta_tip_1/(sgm*Cl_a)))-1);     % λ for each location of t
        lambda_collect_1 =  lambda_collect_1 + lambda_t_1;                      % Collecting all λ over one gaussian span for mean
        dCT_dr_1(j) = 0.5*sgm*Cl_a*(theta_tip_1 - lambda_t_1)*t(j);             % CT distribution from lambda (Q1-d)
        f_ct_1(j) = dCT_dr_1(j);                                                % Gaussian integration function for CT
        CT_1(i) = CT_1(i) + W(j)*f_ct_1(j)*((b-a)/2);                           % Numerically integration of CT over t
        f_cp_1(j) = dCT_dr_1(j)*lambda_t_1;                                     % Gaussian integration function for CPi
        % f_cp_1(j) = 0.5*sgm*Cl_a*(theta_tip_1 - lambda_t_1)*t(j)*lambda_t_1;
        CPi_1(i) = CPi_1(i) + W(j)*f_cp_1(j)*((b-a)/2);                         % Numerically integration of CPi over t
        f_cp_p1(j) = 0.5*sgm*Cd_0*(t(j)^3);                                     % Gaussian integration function for CPp
        CPp_1(i) = CPp_1(i) + W(j)*f_cp_p1(j)*((b-a)/2);

        % Case-2: No twist: θ_tw = 0 deg

        theta_tip_2 = theta_2(j)*t(j);
        lambda_t_2 = (sgm*Cl_a/16)*(sqrt(1+(32*theta_tip_1/sgm/Cl_a))-1);
        lambda_collect_2=  lambda_collect_2 + lambda_t_2;
        dCT_dr_2(j)= 0.5*sgm*Cl_a*(theta_tip_2-lambda_t_2)*t(j);
        f_ct_2(j)= dCT_dr_2(j);
        CT_2(i)= CT_2(i) + W(j)*f_ct_2(j)*0.5*(b-a);
        f_cp_2(j)= 0.5*sgm*Cl_a*(theta_tip_2-lambda_t_2)*t(j)*lambda_t_2;
        CPi_2(i)= CPi_2(i) + W(j)*f_cp_2(j)*0.5*(b-a);
        f_cp_p2(j)=0.5*sgm*Cd_0*(t(j)^3);
        CPp_2(i)=CPp_2(i)+W(j)*f_cp_p2(j)*0.5*(b-a);
        
        % Case-3: Linear Twist: θ_tw = -15 deg

        theta_tip_3 = theta_3(j)*t(j);
        lambda_t_3 = (sgm*Cl_a/16)*((1+(32*theta_tip_3/sgm/Cl_a))^0.5-1);
        lambda_collect_3 =  lambda_collect_3 + lambda_t_3;
        dCT_dr_3(j) = 0.5*sgm*Cl_a*(theta_tip_3-lambda_t_3)*t(j);
        f_ct_3(j) = dCT_dr_3(j);
        CT_3(i) = CT_3(i) + W(j)*f_ct_3(j)*0.5*(b-a);
        f_cp_3(j) = 0.5*sgm*Cl_a*(theta_tip_3-lambda_t_3)*t(j)*lambda_t_3;
        CPi_3(i) = CPi_3(i) + W(j)*f_cp_3(j)*0.5*(b-a);
        f_cp_p3(j) = 0.5*sgm*Cd_0*(t(j)^3);
        CPp_3(i) = CPp_3(i)+W(j)*f_cp_p3(j)*0.5*(b-a);
        
        % Case-4: Ideal Twist: θ_tip and Ideal Taper: c_tip = 0.5*c

        theta_tip_4 = theta_tip;
        lambda_t_4 = (sgm_4(j)*Cl_a/16)*((1+(32*theta_tip_4/sgm_4(j)/Cl_a))^0.5-1);
        lambda_collect_4 =  lambda_collect_4 + lambda_t_4;
        dCT_dr_4(j) = 0.5*sgm_4(j)*Cl_a*(theta_tip_4-lambda_t_4)*t(j);
        f_ct_4(j) = dCT_dr_4(j);
        CT_4(i) = CT_4(i) + W(j)*f_ct_4(j)*0.5*(b-a);
        f_cp_4(j) = 0.5*sgm_4(j)*Cl_a*(theta_tip_4-lambda_t_4)*t(j)*lambda_t_4;
        CPi_4(i) = CPi_4(i) + W(j)*f_cp_4(j)*0.5*(b-a);
        f_cp_p4(j) =0.5*sgm_4(j)*Cd_0*(t(j)^3);
        CPp_4(i) = CPp_4(i)+W(j)*f_cp_p4(j)*0.5*(b-a);
    end
    
    loc(k)= t(3); % for each gaussian span it is a middle location: for plotting the mean values over the entire span

    % Mean values for each collected Non-dimentional parameter in Case-1

    lambda_mean_1(k) = lambda_collect_1/6;                  % Induced inflow
    lambda_collect_1 = 0;
    CT_1_mean(k) = mean(dCT_dr_1);                          % Thrust distribution    
    CPi_1_mean(k) = CPi_1(i);                               % Induced power
    CPp_1_mean(k) = CPp_1(i);                               % Profile power
    CP_1_mean(k) = CPi_1_mean(k) + CPp_1_mean(k);           % Total power

    % Mean values for each collected Non-dimentional parameter in Case-2

    lambda_mean_2(k)= lambda_collect_2/6;
    lambda_collect_2=0;
    CT_2_mean(k)= mean(dCT_dr_2);
    CPi_2_mean(k)=CPi_2(i);
    CPp_2_mean(k)=CPp_2(i);
    CP_2_mean(k)= CPi_2_mean(k) + CPp_2_mean(k);

    % mean values for each collected Non-dimentional parameter in Case-3

    lambda_mean_3(k)= lambda_collect_3/6;
    lambda_collect_3=0;
    CT_3_mean(k)= mean(dCT_dr_3);
    CPi_3_mean(k)=CPi_3(i);
    CPp_3_mean(k)=CPp_3(i);
    CP_3_mean(k)= CPi_3_mean(k) + CPp_3_mean(k);

    % mean values for each collected Non-dimentional parameter in Case-4

    lambda_mean_4(k)= lambda_collect_4/6;
    lambda_collect_4=0;
    CT_4_mean(k)= mean(dCT_dr_4);
    CPi_4_mean(k)=CPi_4(i);
    CPp_4_mean(k)=CPp_4(i);
    CP_4_mean(k)= CPi_4_mean(k) + CPp_4_mean(k);
    
    k = k+1; % Incrementing the location
    
end

results_BEMT_Q1
plot_BEMT_Q1

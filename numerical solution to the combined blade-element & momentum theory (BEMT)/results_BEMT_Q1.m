%% Results

r = R_cut:(R_cut/10):1;

% Ideal Twist
CT_1_total = sum(CT_1);                                                     % Total CT over the entire disc
CP_1_total = sum(CPi_1) + sum(CPp_1);                                       % Total CP over the entire disc
lambda_1 = sqrt(CT_1_total/2);                                              % Recalculating induced inflow from final CT
theta_75_1 = ((6*CT_1_total)/(sgm*Cl_a)) + ((3/2)*lambda_1);                % θ75 for ideal twist
theta_tip_IdealTwist = (4*CT_1_total/(sgm*Cl_a)) + lambda_1;                % θ_tip closed form solution
theta_result_1 = (theta_tip_IdealTwist./r)*(180/pi);                        % θ vs r plot

% No twist: θ_tw = 0 deg
CT_2_total = sum(CT_2);
CP_2_total = sum(CPi_2) + sum(CPp_2);
lambda_2 = sqrt(CT_2_total/2);
theta_result_2 = ((6*CT_2_total/(sgm*Cl_a)) + (3*lambda_2)/2)*(180/pi).*ones(1,length(r));
% theta_result_2= ((6*CT_2_total/sgm/Cl_a)+(3*lambda_2)/2)*(180/pi);
theta_tip_NoTwist = theta_result_2.*r;

% Linear Twist: θ_tw= -15deg
CT_3_total= sum(CT_3);
CP_3_total= sum(CPi_3)+sum(CPp_3);
lambda_3= sqrt(CT_3_total/2);
theta_0_3= ((6*CT_3_total/(sgm*Cl_a))+(3*lambda_3)/2);
theta_result_3= (theta_0_3 + theta_tw_3*(r-0.75))*(180/pi);
theta_tip_Twist= theta_result_3.*r;

% Ideal Twist: θ_tip and Ideal Taper: c_tip= 0.5*c
CT_4_total= sum(CT_4);
CP_4_total= sum(CPi_4)+sum(CPp_4);
lambda_4= (CT_4_total/2)^0.5;
sgm_75= sgm_tip/0.75;
theta_75_4= ((6*CT_4_total/(sgm_75*Cl_a))+(3*lambda_4)/2);
theta_result_4= (((4*CT_4_total/(sgm_tip*Cl_a))+lambda_4)./r)*(180/pi);
theta_tip_IdealTaper= (4*CT_4_total/(sgm_tip*Cl_a))+lambda_4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


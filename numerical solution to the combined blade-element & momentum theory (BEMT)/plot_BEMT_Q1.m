%% ## Plot 1  'Pitch Angle v/s r'

figure('Name','Pitch Angle v/s r','NumberTitle','off')
plot(r,theta_result_1,'r.-',r,theta_result_2,'o-b',r,theta_result_3,'x-g',r,theta_result_4,'k-*')
title({'Variation of Pitch angle (θ) with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Pitch angle (θ)')
xlim([0.2,1])
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper')
grid on

%% ## Plot 2 'Numerical Solution vs. Analytical Solution'

figure('Name','Numerical Solution v/s Analytical Solution','NumberTitle','off')
plot(r,theta_result_1,'k-*',r,(theta_tip./r),'r.-') % Numerical Solution, Analytical formula
title({'Comparison of Numerical solution with','closed form exact formula (for ideal twist)'})
xlabel('Non-dimensional radial location (r)')
ylabel('Pitch angle (θ)')
xlim([0.2,1])
legend('Numerical Solution','Closed form exact formula')
grid on

%% ## Plot 3  Angle of attack v/s non dim. radial location

figure('Name','Variation of AOA v/s r','NumberTitle','off')

% Ideal Twist: θ_tip
alpha_1 = theta_result_1 - (lambda_1./r);      % Calculating AOA for all the r locations

% No twist: θ_tw = 0deg
alpha_2 = theta_result_2 - (lambda_2./r);

% Linear Twist: θ_tw = -15deg
alpha_3 = theta_result_3 - (lambda_3./r);

% Ideal Twist: θ_tip and Ideal Taper: c_tip = 0.5*c
alpha_4 = theta_result_4 - (lambda_4./r);

plot(r,alpha_1,'r.-',r,alpha_2,'o-b',r,alpha_3,'x-g',r,alpha_4,'k-*') 
title({'Variation of Angle of Attack (α) with','non-dimensional radial location'})
% title({'First line','Second line'})
xlabel('Non-dimensional radial location (r)')
ylabel('Angle of Attack (α) [in degrees]')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper')
xlim([0.2,1])
grid on

%% ## Plot 4 Inflow vs r

figure('Name','λ v/s r','NumberTitle','off')
plot(loc,lambda_mean_1,'r.-',loc,lambda_mean_2,'o-b',loc,lambda_mean_3,'x-g',loc,lambda_mean_4,'k-*') % plotting the induced inflow for all r locations 
title({'Variation of Induced Inflow with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Induced Inflow (λ)')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper')
xlim([0.2, 1])
grid on

%% ## Plot 5 dCT/dr vs r

figure('Name','dCT/dr v/s r','NumberTitle','off')
plot(loc,CT_1_mean,'r.-',loc,CT_2_mean,'o-b',loc,CT_3_mean,'x-g',loc,CT_4_mean,'k-*') % plotting the CT distribution for all r locations
title({'Variation of thrust distribution with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Thrust distribution (dC_T/dr)')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper','Location','northwest')
xlim([0.2, 1])
grid on

%% ## Plot 6 CQ vs. r

figure('Name','CQ v/s r','NumberTitle','off')
plot(loc,CP_1_mean,'r.-',loc,CP_2_mean,'o-b',loc,CP_3_mean,'x-g',loc,CP_4_mean,'k-*') % plotting the total power distn. for all r locations
title({'Variation of Total Torque distribution with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Total Torque distribution (C_Q)')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper','Location','northwest')
xlim([0.2 1])
grid on

%% ## Plot 7 CQi v/s r

figure('Name','CQi v/s r','NumberTitle','off')
plot(loc,CPi_1_mean,'r.-',loc,CPi_2_mean,'o-b',loc,CPi_3_mean,'x-g',loc,CPi_4_mean,'k-*') % plotting the induced power distn. for all r locations
title({'Variation of Induced Torque distribution with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Induced torque distribution (C_Q_i)')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper','Location','northwest')
xlim([0.2, 1])
grid on

%% ## Plot 8 CQp vs r

figure('Name','CQp v/s r','NumberTitle','off')
plot(loc,CPp_1_mean,'r.-',loc,CPp_2_mean,'o-b',loc,CPp_3_mean,'x-g',loc,CPp_4_mean,'k-*') % plotting the profile power distn. for all r locations
title({'Variation of Profile Torque distribution with','non-dimensional radial location'})
xlabel('Non-dimensional radial location (r)')
ylabel('Profile Torque distribution (C_Q_p)')
legend('Ideal Twist : θ_t_i_p','No twist : θ_t_w = 0 deg'...
    ,'Linear Twist : θ_t_w = -15 deg','Ideal Twist : θ_t_i_p & Ideal Taper','Location','northwest')
xlim([0.2,1])
grid on

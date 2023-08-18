%% Calculation of inflow ratio at the rotor disk in forward flight from the momentum theory and Solution of λ using Fixed Point Iteration compared with the approximate (valid for µ > 0.2) solution of λ

clear all
close all
clc

% Given Data
alpha = 2*pi/180;               % Disk Angle of Attack (α)
Ct = 0.006;                     % Thrust coefficient

%% Newton Raphson
lambda_0 = sqrt(Ct/2);          % Initial guess of inflow ratio (typically hover inflow ratio)
lambda(1) = lambda_0;

for i = 1:1:301
    mu(i) = 0.002*i - 0.002;                                                % 0<=µ<=0.6 ~~ varying Advanced ratio (Incompressible Subsonic regime)
    lambda_apprx(i) = mu(i)*tan(alpha) + Ct/(2*mu(i));
    f = lambda(i) - mu(i)*tan(alpha) - 0.5*Ct/sqrt(mu(i)^2 + lambda(i)^2);
    df = 1 + 0.5*(Ct*lambda(i))/((mu(i)^2 + lambda(i)^2)^1.5);
    lambda(i+1) = lambda(i) - (f/df);
    lambda_nr(i) = lambda(i+1);
end

% Plotting of data (Newton Raphson)
plot(mu,lambda_nr,'r')
hold on
plot(mu,lambda_apprx,'--k')
xlim([0 0.6])
ylim([0 0.08])
title('µ v/s λ (Newton Raphson method & Approximate Solution)')
xlabel('Advance Ratio, µ')
ylabel('Inflow Ratio, λ')
legend('Newton Raphson','Approximate Solution')
grid on
hold off
annotation('textbox', [0.55, 0.68, 0.34, 0.1], 'String', "Disk Angle of Attack (α) = 2°, Thrust coefficient (Ct) = 0.006")

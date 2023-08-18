%% Calculation of the inflow ratio at the rotor disk in forward flight from momentum theory

clear all
close all
clc

% Given Data
alpha = 2*pi/180;               % Disk Angle of Attack (α)
mu = 0.2;                       % Advanced ratio (µ)
Ct = 0.006;                     % Thrust coefficient

%% Fixed Point Iteration
lambda_0 = sqrt(Ct/2);          % Initial guess of inflow ratio (typically hover inflow ratio)
lambda = lambda_0;
lambda_old = lambda_0;          % Old value of Inflow ratio
maxerrper = 0.0001;             % Maximum error
n = 10;                         % No. of iterations

for i=1:n
    lambda = mu*tan(alpha) + Ct/(2*sqrt(mu^2 + lambda^2));
    errper(i) = abs((lambda - lambda_old)/lambda_old)*100;
    errper_fpi(i) = errper(i);
    if (errper_fpi(i) < maxerrper)
        break;
    end
    lambda_old = lambda;
    i = i+1;
end

% Plotting of data (Fixed Point Iteration)
plot(1:i,errper_fpi,'bo',1:i,errper_fpi,'-r')
xlim([0,6])
ylim([0,70])
title('Inflow ratio (λ) calculation using Momentum theory')
ylabel('% error');
xlabel('No of iteration');
grid on
hold on

%% Newton Raphson Iteration
lambda_0 = sqrt(Ct/2);          % Initial guess of inflow ratio (typically hover inflow ratio)
lambda = lambda_0;
lambda_old = lambda_0;          % Old value of Inflow ratio
maxerrper = 0.0001;             % Maximum error
n = 10;                         % No. of iterations

for i=1:n
    f = lambda - mu*tan(alpha) - 0.5*Ct/sqrt(mu^2 + lambda^2);
    df = 1 + 0.5*(Ct*lambda)/((mu^2 + lambda^2)^1.5);
    lambda = lambda - (f/df);
    errper(i) = abs((lambda - lambda_old)/lambda_old)*100;
    errper_nr(i) = errper(i);
    if (errper_nr(i) < maxerrper)
        break;
    end
    lambda_old = lambda;
    i = i+1;
end

% Plotting of data (Newton Raphson Iteration)
plot(1:i,errper_nr,'k*',1:i,errper_nr,'--k')
legend({'λ (Fixed Point Iteration)','','λ (Newton Raphson)'},'Location','northeast')
hold off
annotation('textbox', [0.55, 0.65, 0.35, 0.15], 'String', "Disk Angle of Attack (α) = 2°, Advanced ratio (µ) = 0.2, Thrust coefficient (Ct) = 0.006")


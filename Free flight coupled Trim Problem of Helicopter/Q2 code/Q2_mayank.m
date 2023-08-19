%%Coupled trim solution for the helicopter
clear all
close all;
clc
%%Given Input Data:
Nb=4;
rho =0.002377*32.174;       % unit in FPS
R=26.83;                    % unit in FPS
c=1.509;                    % unit in FPS
omega =27;
v_tip=omega*R;
sigma=Nb*c/(pi*R);
step=0.02;
mu(1)=0;
i=1;
while mu(i)<=0.42
    [theta_o(i), theta_1c(i), theta_1s(i),alpha_s(i),phi_s(i),lamda(i)]=Controlinput(mu(i)); 
    fndms = Q1func(theta_o(i), theta_1c(i), theta_1s(i), lamda(i),mu(i), phi_s(i), alpha_s(i));
    Cmx_sig(i)=fndms(4)/(rho*pi*(R^2)*(v_tip^2)*R*sigma);
    Cmy_sig(i)=fndms(5)/(rho*pi*(R^2)*(v_tip^2)*R*sigma);
    mu(i+1)= mu(i)+step
    i=i+1
end
%%
mu(end)=[];
% set(0,'DefaultAxesFontSize',16);
% set(0,'DefaultAxesFontWeight','bold');
% set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on');
figure(1)
plot(mu,theta_o*180/pi,'r')
hold on
plot(mu,theta_1c*180/pi,'b')
hold on
plot(mu,theta_1s*180/pi,'g')
hold off
xlabel('\mu (Advance ratio)', 'FontWeight', 'bold')
ylabel('Control input angles [degrees]', 'FontWeight', 'bold')
title('Variation of Control input angles v/s \mu', 'FontWeight', 'bold');
legend('\theta_0 (Collective)','\theta_1_c (Lateral cyclic)','\theta_1_s (Longitudinal cyclic)');
xlim([0,0.45]);
saveas(gcf,'fig1.fig');

figure(2)
plot(mu,alpha_s*180/pi,'r')
hold on
plot(mu,phi_s*180/pi,'b')
hold off
xlabel('\mu (Advance ratio)', 'FontWeight', 'bold')
ylabel('Vehicle shaft angles [degrees]', 'FontWeight', 'bold')
title('Variation of Vehicle shaft angles v/s \mu', 'FontWeight', 'bold');
legend('\alpha_s (Longitudinal shaft angle)','\phi_s (Lateral shaft angle)');
xlim([0,0.45]);
saveas(gcf,'fig2.fig');

figure(3)
plot(mu,Cmx_sig)
hold on
plot(mu,Cmy_sig)
xlabel('\mu', 'FontWeight', 'bold')
ylabel('Non dim. mean hub loads', 'FontWeight', 'bold')
title('Non dim. mean hub loads vs \mu');
legend('C_{MX}/\sigma','C_{MY}/\sigma');
xlim([0,0.45]);
%%
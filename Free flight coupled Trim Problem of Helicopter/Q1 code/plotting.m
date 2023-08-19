psi_plot = 1:1:360;

figure(1)
plot(psi_plot,fx(1,:),'r');
hold on
plot(psi_plot,fy(1,:),'g');
hold on
plot(psi_plot,fz(1,:),'b');
hold off
xlabel('\psi [deg]', 'FontWeight', 'bold')
ylabel('Blade hub Shear Forces in rotating frame [N]', 'FontWeight', 'bold')
title({'Blade hub shear forces in the rotating frame'}, {'as a function of azimuth \psi'}, 'FontWeight', 'bold')
legend('fx','fy','fz')
xlim([0 360])
ylim([-50000 300000])
xticks(0:90:360)
xticklabels({'0','90\circ','180\circ','270\circ','360\circ'})
% grid on
% saveas(figure(1),Blade hub shear forces.jpg);
saveas(gcf,'fig1.fig');

figure(2)
plot(psi_plot,mx(1,:),'r');
hold on
plot(psi_plot,my(1,:),'g');
hold on
plot(psi_plot,mz(1,:),'b');
hold off
xlabel('\psi [deg]', 'FontWeight', 'bold', 'FontWeight', 'bold')
ylabel('Blade hub Bending Moments in rotating frame [N-m]', 'FontWeight', 'bold')
title({'Blade hub bending moments in the rotating frame'}, {'as a function of azimuth \psi'}', 'FontWeight', 'bold')
legend('mx','my','mz')
xlim([0 360])
xticks(0:90:360)
xticklabels({'0','90\circ','180\circ','270\circ','360\circ'})
% grid on
saveas(gcf,'fig2.fig');

figure(3)
% Non-dimentional shear forces
fz_nondim_1 = fz(1,:)/(rho*A*(omega*R)^2);
fz_nondim_2 = fz(2,:)/(rho*A*(omega*R)^2);
fz_nondim_3 = fz(3,:)/(rho*A*(omega*R)^2);
fz_nondim_4 = fz(4,:)/(rho*A*(omega*R)^2);
plot(psi_plot,fz_nondim_1,'r')
hold on
plot(psi_plot,fz_nondim_2,'g')
hold on
plot(psi_plot,fz_nondim_3,'b')
hold on
plot(psi_plot,fz_nondim_4,'k')
hold off
xlabel('\psi [deg]', 'FontWeight', 'bold', 'FontWeight', 'bold')
ylabel('Non-dimensional Vertical shear force in rotating frame [N]', 'FontWeight', 'bold')
title({'The variation of non-dimensional rotating frame vertical shear force'}, {'as a function of \psi for all four blades'}, 'FontWeight', 'bold')
legend('Blade 1 (\psi_0=0^o)', 'Blade 2 (\psi_0=90^o)','Blade 3 (\psi_0=180^o)','Blade 4 (\psi_0=270^o)');
xlim([0 360])
xticks(0:90:360)
xticklabels({'0','90\circ','180\circ','270\circ','360\circ'})
% grid on
saveas(gcf,'fig3.fig');

figure(4)
plot(psi_plot,Fx_plot_fixed,'r');
hold on
plot(psi_plot,Fy_plot_fixed,'g');
hold on
plot(psi_plot,Fz_plot_fixed,'b');
hold off
xlabel('\psi [deg]', 'FontWeight', 'bold', 'FontWeight', 'bold')
ylabel('Blade hub Shear Forces in fixed frame [N]', 'FontWeight', 'bold')
title({'Variation of forces in fixed frame of reference'}, {'as a function of \psi'}, 'FontWeight', 'bold')
legend('Fx or H','Fy or Y','Fz or T')
xlim([0 360])
ylim([-10000 90000])
% ylim([-2000 2000])
xticks(0:90:360)
xticklabels({'0','90\circ','180\circ','270\circ','360\circ'})
% grid on
saveas(gcf,'fig4.fig');

figure(5)
plot(psi_plot,Mx_plot_fixed,'r');
hold on
plot(psi_plot,My_plot_fixed,'g');
hold on
plot(psi_plot,Mz_plot_fixed,'b');
hold off
xlabel('\psi [deg]', 'FontWeight', 'bold', 'FontWeight', 'bold')
ylabel('Blade hub Bending Moments in fixed frame [N-m]', 'FontWeight', 'bold')
% ylabel('Y Label', 'FontWeight','bold')
title({'Variation of Moments in fixed frame'}, {'of reference as a function of \psi'}, 'FontWeight', 'bold')
legend('Mx','My','Mz or Q')
xlim([0 360])
% ylim([-2000 2000])
xticks(0:90:360)
xticklabels({'0','90\circ','180\circ','270\circ','360\circ'})
% grid on
saveas(gcf,'fig5.fig');

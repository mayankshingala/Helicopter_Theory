%% Plots

psi = 0 : d_psi : rev*2*pi;

figure('Name','Question-01 (Assignment-04)','NumberTitle','off')
plot(psi,M_beta,'k--.')
hold on
plot(psi,M_beta_exact,'r-')

xlabel('\bf \psi, Azimuth location [rad]')
ylabel('$\bar{\bf M_{\beta}}$, Non-dimensional Aerodynamic flap moment','Interpreter','Latex')

title('Aerodynamic flap moment as a function of blade azimuth');

legend('Numerical Solution', 'Exact Analytical solution')
xlim([0 rev*2*pi])
xticks(0:pi:rev*2*pi)
xticklabels({'\0','\pi','2\pi','3\pi','4\pi','5\pi','6\pi','7\pi','8\pi','9\pi','10\pi'})
grid on
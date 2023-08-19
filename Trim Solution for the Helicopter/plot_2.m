
figure(1)
% figure('Name','Question-01 (Assignment-04)','NumberTitle','off')
xlabel('\bf \psi, Azimuth location [rad]')
ylabel('\bf \beta, Flap angle [rad]')
title('Flap dynamics solution over blade azimuth locations')
xlim([0 2*pi])
xticks(0 : pi/4: 2*pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
grid on
% plotstyle = {'-*','-.','-o','-d','-h','-p'};
plot(psi, beta_HB(1:length(psi),j), 'r-.')
legend({'Time Step: 2\circ', 'Time Step: 3\circ', 'Time Step: 5\circ', ...
        'Time Step: 10\circ', 'Time Step: 15\circ', 'Harmonic Balance'},'Location','north')
hold off

% figure(2)
figure('Name','Question-02-a (Assignment-04)','NumberTitle','off')
psi = 0 : d_psi(1) : rev*2*pi-d_psi(1);
plot(psi, beta(1:length(psi),1), 'r')
hold on
plot(psi, beta_HB(1:length(psi),1), 'b-.')
hold off
xlabel('\bf \psi, Azimuth location [rad]')
ylabel('\bf \beta, Flap angle [rad]')
title({'Comparison of Flap equation solution for One Revolution:', 'Numerical vs Harmonic Balance'})
legend({'Numerical Solution', 'Harmonic Balance Solution'},'Location','north')
xlim([0 2*pi])
xticks(0 : pi/2 : 2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
grid on

% figure(3)
figure('Name','Question-02-b (Assignment-04)','NumberTitle','off')
psi = 0 : d_psi(1) : rev*2*pi-d_psi(1);
plot(psi, beta(1:length(psi),1), 'r')
hold on
plot(psi, beta_HB(1:length(psi),1), 'b-.')
hold off
xlabel('\bf \psi, Azimuth location [rad]')
ylabel('\bf \beta, Flap angle [rad]')
title({'Comparison of Flap equation solution over time','history of 5 revolution: Numerical vs Harmonic Balance'})
ax.TitleHorizontalAlignment = 'left';
legend({'Numerical Solution', 'Harmonic Balance Solution'}, 'Location', 'northeastoutside')
xlim([0 rev*2*pi])
xticks(0 : pi : rev*2*pi)
xticklabels({'0','\pi','2\pi','3\pi','4\pi','7\pi','6\pi','7\pi','8\pi','9\pi','10\pi'})
grid on

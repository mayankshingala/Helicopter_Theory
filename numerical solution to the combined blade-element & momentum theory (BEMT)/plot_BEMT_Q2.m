%% Plot 1: Pitch v/s r
figure('Name','Refrence Pitch Angle (θ) v/s r for No twist case (Hovering flight)','NumberTitle','off')
plot(r,theta75_result,'g.-');
hold on
plot(r,theta75_result_TL,'r-','LineWidth',1.2);
hold off
title('The Reference Pitch input (θ_7_5) needed for hovering flight')
xlabel('Non-dimensional radial location (r)')
ylabel('Pitch angle (θ) [in degrees]')
ylim([0.16 0.19])
legend('No Tip-loss','With Tip-loss (Prandtl tip loss function)','Location','northeast')
grid on

%% plot 2: λ v/s r
figure('Name','λ v/s r for No twist case (Hovering flight)','NumberTitle','off')
plot(loc,lambda_mean,'g.-','LineWidth',1.2)
hold on
plot(loc,lambda_mean_TL,'r-')
hold off
title('The effect of the tip-loss on the induced inflow')
xlabel('Non-dimensional radial location (r)')
ylabel('Induced Inflow (λ)')
xlim([0.2 1])
legend('No Tip-loss','With Tip-loss (Prandtl tip loss function)','Location','northwest')
grid on

%% Plot 3: CT distribution v/s r 
figure('Name','dCT/dr v/s r for No twist case (Hovering flight)','NumberTitle','off')
plot(loc,CT_mean,'g.-','LineWidth',1.2)
hold on
plot(loc,CT_TL_mean,'r-')
hold off
title('The effect of the tip-loss on the Thrust distribution')
xlabel('Non-dimensional radial location (r)')
ylabel('Thrust distribution (C_T)')
legend('No Tip-loss','With Tip-loss (Prandtl tip loss function)','Location','northwest')
grid on

%% Plot 4: CP distribution v/s r
figure('Name','dCQ/dr v/s r for No twist case (Hovering flight)','NumberTitle','off')
plot(loc,CP_mean,'g.-','LineWidth',1.2)
hold on
plot(loc,CP_TL_mean,'r-')
hold off
title('The effect of the tip-loss on the Torque distribution')
xlabel('Non-dimensional radial location (r)')
ylabel('Total Torque distribution (C_Q)')
legend('No Tip-loss','With Tip-loss (Prandtl tip loss function)','Location','northwest')
grid on

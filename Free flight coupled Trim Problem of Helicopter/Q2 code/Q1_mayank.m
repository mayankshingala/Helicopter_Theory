%%Coupled trim solution for the helicopter
clc
clear all
close all
rho = 1.225;        %kg/m3 Density of air
Nb = 4;             %Number of blades
R = 8.18;           %m  Blade radius
c = 0.46;           %m  Blade chord
Cdo = 0.01;         %Profile Drag coefficient
Cl_a = 5.73;        %Lift curve slope
omega = 27 ;        % rad/sec  Rotor angular speed
gamma_b = 1.04;     %per rev  Blade flap frequency
gamma = 8.0;        %Lock number,
W = 70000 ;         %N  Weight of the aircraft
theta_tw = 0*pi/180;%Blade twist rate,
g = 9.81;
e = 0;              %non dimnesional root cutout
Cd = Cdo;           % assumption
A = pi*R^2;

%input values
mu = 0.2857; lambda = 0.02284;
theta_0 = 8.2*pi/180; theta_1c = 3.3*pi/180; theta_1s = -11.2*pi/180;
alpha_s = -4.1*pi/180;  phi_s= -3.84*pi/180;
%%
step = 1;
psir = (0:step:360)*pi/180;
dpsi = psir(2) - psir(1);
psiall = [psir; psir+(pi/2); psir+(pi); psir+(3*pi/2)];

%% Harmonic balance
beta_0 = (gamma/(gamma_b^2))*((theta_0*(1+(mu^2))/8) + (theta_tw*(1 + (5*(mu^2)/6))/10) + (mu*theta_1s/6) - (lambda/6));
consta = (gamma_b^2) - 1;
constb = gamma*((theta_1c*(1+(0.5*(mu^2)))/8) - (mu*beta_0/6));
constc = gamma*((1+(0.5*(mu^2)))/8);
constd = gamma*((theta_1s*(1-(0.5*(mu^2)))/8) + (mu*theta_0/3) - (mu*lambda/4) + ((mu^2)*theta_1s/4) + (mu*theta_tw/4));
conste = -(1/8)*(gamma)*(1-(0.5*(mu^2)));
matA = [consta, constc; conste, consta];
matB = [constb; constd];
[beta_const]  = inv(matA)*matB;
beta_1c = beta_const(1);
beta_1s = beta_const(2);
num=1;
for psi = (0:step:360)*pi/180
    
    beta(num) = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);
    beta_d(num) = (-beta_1c*sin(psi) + beta_1s*cos(psi));
    beta_dd(num) = 0;
    
    num = num+1;
    
end



Ct = W/(rho*A*(omega*R)^2);
lambda_0 = sqrt(Ct/2);
lambda = sqrt(Ct/2);
error =1;

while error>0.01
    
    lambda = (mu*tan(alpha_s)) + (Ct/(2*sqrt(mu^2 + lambda_0^2)));
    error = lambda-lambda_0;
    lambda_0 = lambda;
    
    
end
%% Guassian Quadrature and NewMarks Method to solve flap eqn
x(1)=-0.932469514; x(6)=-x(1);
x(2)=-0.661209386; x(5)=-x(2);
x(3)=-0.0238619186; x(4)=-x(3);
w(1)=0.171324492; w(6)=w(1);
w(2)=0.360761573; w(5)=w(2);
w(3)=0.467913935; w(4)=w(3);
a = linspace(0.15,1,11);
[r, c] = size(psiall);
fx = zeros(r,c-1);
fy = zeros(r,c-1);
fz = zeros(r,c-1);

mx = zeros(r,c-1);
my = zeros(r,c-1);
mz = zeros(r,c-1);

Fx_nrot = zeros(r,c-1);
Fy_nrot = zeros(r,c-1);
Fz_nrot = zeros(r,c-1);

Mx_nrot = zeros(r,c-1);
My_nrot = zeros(r,c-1);
Mz_nrot = zeros(r,c-1);



for all=1:1
    psi = psiall(all,:);
    
    
    for num = 1:length(psi)-1
        
        for i=1:10
            
            for j=1:6
                
                t(j) = ( (a(i+1)-a(i))*x(j)/2 ) + ( (a(i)+a(i+1))/2 );
                
                ut(j) = t(j) + mu*sin(psi(num));
                up(j) = lambda*cos(beta(num)) + t(j)*beta_d(num) + mu*cos(psi(num))*sin(beta(num))*omega*R;
                
                theta(num) = theta_0 + theta_1c*cos(psi(num)) + theta_1s*sin(psi(num)) + theta_tw*x(j);
                
                
                %for forces (dimensional)
                
                Ut(j) = omega*t(j) + mu*sin(psi(num))*omega*R;
                Up(j) = (lambda*cos(beta(num))*omega*R) + (t(j)*beta_d(num)) + (mu*cos(psi(num))*sin(beta(num))*omega*R);
                Ur(j) = mu*omega*R*cos(psi(num));
                
                
                dFz(j) = 0.5*rho*c*Cl_a*(((Ut(j)^2)*theta(num)) - (Up(j)*Ut(j)))*w(j)*((a(j+1)-a(j))/2);
                dFx(j) = 0.5*rho*c*Cl_a*((Up(j)*Ut(j)*theta(num)) - (Up(j)^2) + (Cd*(Ut(j)^2)/Cl_a))*w(j)*((a(j+1)-a(j))/2);
                dFr(j) = -0.5*beta(num)*Cl_a*rho*c*(((Ut(j)^2)*theta(num)) - (Up(j)*Ut(j)))*w(j)*((a(j+1)-a(j))/2);
                
                dFz_dr(j) = 0.5*rho*c*Cl_a*(((Ut(j)^2)*theta(num)) - (Up(j)*Ut(j)));
                dFx_dr(j) = 0.5*rho*c*Cl_a*((Up(j)*Ut(j)*theta(num)) - (Up(j)^2) + (Cd*(Ut(j)^2)/Cl_a));
                dFr_dr(j) = -0.5*beta(num)*Cl_a*rho*c*(((Ut(j)^2)*theta(num)) - (Up(j)*Ut(j)));
                
                % dimensinal blade shear force
                Sz_rot(j) = (dFz_dr(j))*w(j)*((a(j+1)-a(j))/2);
                Sx_rot(j) = dFx_dr(j)*w(j)*((a(j+1)-a(j))/2);
                Sr_rot(j) = ( - (beta(num)*dFr_dr(j)))*w(j)*((a(j+1)-a(j))/2);
                % dimensional blade moments
                mf_rot(j) = ((t(j)-e)*dFz_dr(j))*w(j)*((a(j+1)-a(j))/2);
                ml_rot(j) = ((t(j)-e)*dFx_dr(j))*w(j)*((a(j+1)-a(j))/2);
                mt_rot(j) = 0;             % assuming zero pitching moment
            end
            
            Sx_rot_sec(i) = sum(Sx_rot);
            Sy_rot_sec(i) = sum(Sr_rot);
            Sz_rot_sec(i) = sum(Sz_rot);    
            mf_rot_sec(i) = sum(mf_rot);
            mt_rot_sec(i) = sum(mt_rot);
            ml_rot_sec(i) = sum(ml_rot);        
        end  
        Sx(num) = sum(Sx_rot_sec);             %rotating frame blade load for different psi(num)
        Sy(num) = sum(Sy_rot_sec);
        Sz(num) = sum(Sz_rot_sec);
        fx(all,num) = Sx(num);             %rotating frame hub load
        fy(all,num) = Sy(num);
        fz(all,num) = Sz(num);
        mf(num) = sum(mf_rot_sec);     %rotating frame blade moment
        mt(num) = sum(mt_rot_sec);
        ml(num) = sum(ml_rot_sec);
        mx(all,num) = mf(num) + (e*Sz(num));            % dimensional hub moments rotating
        my(all,num) = mt(num);
        mz(all,num) = - ml(num) - (e*Sx(num));
        
    end
end
fx(2,1:270) = fx(1,91:360);
fx(2,271:360) = fx(1,1:90);
fx(3,1:180) = fx(1,181:360);
fx(3,181:360) = fx(1,1:180);
fx(4,1:90) = fx(1,271:360);
fx(4,91:360) = fx(1,1:270);

fy(2,1:270) = fy(1,91:360);
fy(2,271:360) = fy(1,1:90);
fy(3,1:180) = fy(1,181:360);
fy(3,181:360) = fy(1,1:180);
fy(4,1:90) = fy(1,271:360);
fy(4,91:360) = fy(1,1:270);

fz(2,1:270) = fz(1,91:360);
fz(2,271:360) = fz(1,1:90);
fz(3,1:180) = fz(1,181:360);
fz(3,181:360) = fz(1,1:180);
fz(4,1:90) = fz(1,271:360);
fz(4,91:360) = fz(1,1:270);

mx(2,1:270) = mx(1,91:360);
mx(2,271:360) = mx(1,1:90);
mx(3,1:180) = mx(1,181:360);
mx(3,181:360) = mx(1,1:180);
mx(4,1:90) = mx(1,271:360);
mx(4,91:360) = mx(1,1:270);

my(2,1:270) = my(1,91:360);
my(2,271:360) = my(1,1:90);
my(3,1:180) = my(1,181:360);
my(3,181:360) = my(1,1:180);
my(4,1:90) = my(1,271:360);
my(4,91:360) = my(1,1:270);

mz(2,1:270) = mz(1,91:360);
mz(2,271:360) = mz(1,1:90);
mz(3,1:180) = mz(1,181:360);
mz(3,181:360) = mz(1,1:180);
mz(4,1:90) = mz(1,271:360);
mz(4,91:360) = mz(1,1:270);

for all = 1:4
    for num = 1:length(psi)-1   
        % fixed hub loads and moments
        psi_m(all) =  psi(num) + ((all - 1)*2*pi/Nb);  
        Fx_nrot(all,num) = (fy(all,num)*cos(psi_m(all)) + fx(all,num)*sin(psi_m(all)));
        Fy_nrot(all,num) = (fy(all,num)*sin(psi_m(all)) - fx(all,num)*cos(psi_m(all)));
        Fz_nrot(all,num) = (fz(all,num));    
        Mx_nrot(all,num) = (mx(all,num)*sin(psi_m(all)) + my(all,num)*cos(psi_m(all)));
        My_nrot(all,num) = (-mx(all,num)*cos(psi_m(all)) + my(all,num)*sin(psi_m(all)));
        Mz_nrot(all,num) = (mz(all,num));    
    end
    
end
for num = 1:length(psi)-1
    
    % HUb load rotating frame
    ffx_hrot(num) = sum(fx(:,num));    
    ffy_hrot(num) = sum(fy(:,num));  
    ffz_hrot(num) = sum(fz(:,num));

    %Hub moment rotating frame
    mmx_hrot(num) = sum(mx(:,num));
    mmy_hrot(num) = sum(my(:,num));
    mmz_hrot(num) = sum(mz(:,num));
    Fx(num) = sum(Fx_nrot(:,num));          % fixed frame hub load
    Fy(num) = sum(Fy_nrot(:,num));
    Fz(num) = sum(Fz_nrot(:,num)); 
    Mx(num) = sum(Mx_nrot(:,num));          % fixed frame hub moment
    My(num) = sum(My_nrot(:,num));
    Mz(num) = sum(Mz_nrot(:,num));
    
end
psii = (1:step:360);
figure(1)
plot(psii,fx(1,:),psii,fy(1,:),psii,fz(1,:));
xlabel('\psi (in deg)')
ylabel('Blade hub Shear Forces in rotating frame(N)')
title('Blades hub shear loads (Rotating frame)')
legend('X direction component: fx', 'Y direction component: fy', 'Z direction component: fz');
set(gca,'FontSize',14);

figure(2)
plot(psii,mx(2,:),psii,my(2,:),psii,mz(2,:),'-b');
xlabel('\psi (in deg)')
ylabel('Blades hub bending moment in rotating frame (N-m)')
title('Blades hub bending moment (Rotating frame)')
legend('X direction component: mx', 'Y direction component: my', 'Z direction component: mz');
set(gca,'FontSize',14);

figure(3)
plot(psii,fz/(rho*pi*R^2*(omega*R)^2))
hold on
plot(psii,ffz_hrot/(rho*pi*R^2*(omega*R)^2))
xlabel('\psi (in deg)')
ylabel('non-dimensional rotating frame vertical shear force')
legend('Blade 1(psi_0=0^o)', 'Blade 2(psi_0=90^o)','Blade 3(psi_0=180^o)','Blade 4(psi_0=270^o)');
set(gca,'FontSize',14);

figure(4)
plot(psii,Fx,psii,Fy,psii,Fz);
xlabel('\psi (in deg)')
ylabel('Blade hub Shear forces in fixed frame(N)')
title('Hub forces (Fixed frame)')
legend('X direction component: Fx', 'Y direction component: Fy', 'Z direction component: Fz');
set(gca,'FontSize',14);

figure(5)
plot(psii,Mx,psii,My,psii,Mz)
xlabel('\psi (in deg)')
ylabel('Blade hub Moments in fixed frame(N-m)')
title('HUB Moments in fixed frame');
legend('X direction component: Mx', 'Y direction component: My', 'Z direction component: Mz');
set(gca,'FontSize',14);



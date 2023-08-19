
function Res = Residuefunc( theta_0, theta_1c, theta_1s, lamda,mu, phi_s, alpha_s)

W = 15736.626*32.147;
r =26.83727;
xcg = -1.9;
ycg = 0;
h = 6.023;
f = 19.9;
rho =   0.002377*32.174;
M_xf=0;
M_yf=0;
omega=27;
v_tip=omega*r;

fndm = Q1func(theta_0, theta_1c, theta_1s, lamda,mu, phi_s, alpha_s);

T=fndm(3);%sum(Fz)/length(Fz);
H=fndm(1);%sum(Fx)/length(Fx);
Y=fndm(2);%sum(Fy)/length(Fy);

Q=fndm(6);%sum(Mz)/length(Mz);
MX=fndm(4);%sum(Mx)/length(Mx);
MY=fndm(5);%sum(My)/length(My);

D=0.5*rho*(mu*v_tip)^2*f;
Yf=0;
theta_fp=0;
Ct=T/(rho*(pi*r^2)*v_tip^2);

R1=W-T*cos(alpha_s)*cos(phi_s)+Y*sin(phi_s)-H*sin(alpha_s)*cos(phi_s)+Yf*sin(phi_s)+D*sin(theta_fp);
R2=D*cos(theta_fp)+H*cos(alpha_s)-T*sin(alpha_s);
R3=(Y+Yf)*cos(phi_s)+T*cos(alpha_s)*sin(phi_s)+H*sin(alpha_s)*sin(phi_s);
R4=MY+M_yf-W*(xcg*cos(alpha_s)-h*sin(alpha_s))-D*(xcg*sin(alpha_s)+h*cos(alpha_s));
R5=MX+M_xf+Yf*h+W*(h*sin(phi_s)-ycg*cos(phi_s));
R6=lamda-lamdacal(mu,alpha_s,Ct) ; %mu*tan(alpha_s)+Ct/(2*sqrt(mu^2+lambda^2));

Res=[R1;R2;R3;R4;R5;R6];

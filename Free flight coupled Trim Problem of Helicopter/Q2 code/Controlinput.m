
function [theta_0, theta_1c, theta_1s,alpha_s,phi_s,lambda]=Controlinput(mu)
%%%%initial guess%%%%%%
theta_0= 10*pi/180 ;
theta_1c= 1.5*pi/180;
theta_1s = -8*pi/180;
alpha_s=-5*pi/180;
phi_s=-3*pi/180;
lambda=sqrt(0.0066/2);
if(mu==0)                   %hover case
    lambda=sqrt(0.0066/2);
end
if(mu~=0)                    %cases other than hover
    lambda=mu*tan(alpha_s)+0.0066*1.15/(2*sqrt(mu^2+lambda^2));
end
R=Residuefunc(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
error =norm(R);
delta_theta=zeros(6,1);
i=1;
while error>1e-2
    %update control inputs
    theta_0=theta_0+ delta_theta(1);
    theta_1c= theta_1c+ delta_theta(2);
    theta_1s= theta_1s +delta_theta(3) ;
    lambda=lambda+delta_theta(4);    %update lambda
    phi_s=phi_s+delta_theta(5);       %update phi_s
    alpha_s=alpha_s+delta_theta(6);    %update alpha_s
    R=Residuefunc(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s); %compute residual
    J=jacob(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);  %compute Jacobian
    delta_theta=-J\R;
    i=i+1;
    error=norm(R);
end

end





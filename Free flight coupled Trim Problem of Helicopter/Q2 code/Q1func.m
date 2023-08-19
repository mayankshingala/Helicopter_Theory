%%[Force and moment calculation]

function output  = Q1func( theta_0, theta_1c, theta_1s,lambda, mu , phi_s ,alpha_s)
Nb = 4;
R = 26.83 ;
c = 1.5;
cl_alpha =5.73;
Cdo = 0.01;
omega=27;
Vtip=omega*R;
Nu_Beta = 1.05;
lockn= 8;
rho =0.002377*32.174;
theta_tw = 0*pi/180;
M_m = 1;
C_m = 0;
K_m = Nu_Beta^2;
beta_m = 1/4;
gamma_m = 1/8;
beta_n_t = 0;
beta_n1_t = 0;
beta_n2_t = 0;
x(1)=-0.932469514;
x(2)=-0.661209386;
x(3)=-0.0238619186;
x(4)=0.0238619186;
x(5)=0.661209386;
x(6)=0.932469514;
w(1)=0.171324492;
w(2)=0.360761573;
w(3)=0.467913935;
w(4)=0.171324492;
w(5)=0.360761573;
w(6)=0.467913935;

error =0.1;
azm = (0:1:360)*pi/180 ;
delta =  azm(2) -azm(1);
beta_nl = zeros(1,length(azm));
it = 0;
while((error > 0.000001))
    it = it + 1;
    for i=1:length(azm)-1;
        beta_n(i) = beta_n_t;
        beta_n_1(i) = beta_n1_t;
        beta_n_2(i) = beta_n2_t;
        
        %Gaussian quadrature for M_beta
        r = 0:0.1:1;
        for j=1:length(r)-1
            Sum=0;
            for k=1:6
                t(k) = (((r(j+1) - r(j))/2)* x(k) + (r(j+1)+r(j))/2);
                u_t = (t(k) + mu*sin(azm(i)));
                theta = theta_0 + theta_1c*cos(azm(i)) + theta_1s * sin(azm(i)) + theta_tw * t(k);
                u_p = (lambda + t(k)*beta_n_1(i) + mu*beta_n(i)*cos(azm(i)));
                New_Func(k) = 0.5*t(k)*((u_t^2 * theta) - (u_p * u_t));
                Coeff = w(k)*New_Func(k);
                Sum = Sum + Coeff;
            end
            Mb(j) = ((r(j+1)-r(j))/2) * Sum;
        end
        MB(i) = sum(Mb);
        
        %Newmark's Method
        m = M_m + (C_m*delta*gamma_m) + (K_m*(delta^2)*beta_m);
        F = lockn * MB(i);
        B = C_m*(beta_n_1(i) + delta*(1 - gamma_m)*beta_n_2(i)) + K_m*((beta_n(i)+ delta*beta_n_1(i) + ((delta^2/2)*(1 - 2*beta_m)*beta_n_2(i))));
        beta_n_2(i+1) = (F - B)/m;
        beta_n_1(i+1) = beta_n_1(i) + delta*((1-gamma_m)*beta_n_2(i) + (gamma_m*beta_n_2(i+1)));
        beta_n(i+1) = beta_n(i) + delta*beta_n_1(i) + (delta)^2/2*((1-2*beta_m)*beta_n_2(i) + (2*beta_m*beta_n_2(i+1)));
        
        beta_n_t = beta_n(i+1);
        beta_n1_t = beta_n_1(i+1);
        beta_n2_t = beta_n_2(i+1);
        c_b((it-1)*180+i)=beta_n(i+1);
        c_azm((it-1)*180+i)=azm(i)*180/pi+(it-1)*360;
    end
    error = norm(beta_n - beta_nl)/norm(beta_n);
    beta_nl = beta_n;
end

%Harmonic Balance
beta_h_b = zeros(1,length(azm));
beta_o = lockn*(theta_0*(1+mu^2)/8 + theta_tw*(1+5/6*mu^2)/10 + mu*theta_1s/6 - lambda/6)/Nu_Beta^2;
m = [(Nu_Beta^2-1) (1+mu^2/2)*lockn/8 ; -(1-mu^2/2)*lockn/8 (Nu_Beta^2-1)];
n = [lockn*(theta_1c*(1 + mu^2/2)/8 - mu*beta_o/6);...
    lockn*(theta_1s*(1 - mu^2/2)/8 + mu*theta_0/3 - mu*lambda/4 + mu^2*theta_1s/4 + mu*theta_tw/4)];
p = m\n;
beta_1c = p(1);
beta_1s = p(2);
for i = 1:length(azm)
    beta_h_b(i) = beta_o + beta_1c*cos(azm(i)) + beta_1s*sin(azm(i));
    
end
azm = (0:1:360)*pi/180 ;
r = 1:1:25;
beta = beta_n;
beta_dot = beta_n_1;
for i=1:length(azm)-1;
    for j=1:length(r)
        Sum=[0 0 0 0 0 0 0 0 0];
        for k=1:6
            t(k) = 0.5* x(k) + (r(j)+0.5);
            u_t(i) = omega*t(k) + mu*omega*R*sin(azm(i));
            theta(i) = theta_0 + theta_1c*cos(azm(i)) + theta_1s * sin(azm(i)) + theta_tw * t(k);
            u_p(i) = (lambda*omega*R*cos(beta(i)) + omega* t(k)*beta_dot(i) + mu*omega*R*cos(azm(i))*sin(beta(i)));
            fz(k) = w(k)*(0.5)*0.5*rho*c*cl_alpha*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            Sum(1) = Sum(1) + fz(k);
            fx(k) = w(k)*(0.5)*0.5*rho*c*cl_alpha*(((Cdo/cl_alpha)*u_t(i)^2  + (u_p(i) * u_t(i)*theta(i)) -u_p(i)^2 ));
            Sum(2) = Sum(2) + fx(k);
            fr(k) = -w(k)*(0.5)*0.5*rho*c*cl_alpha*beta(i)*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            Sum(3) = Sum(3) + fr(k);
            dfz(k) = 0.5*rho*c*cl_alpha*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            mb(k) =  w(k)*(0.5*t(k)*dfz(k))/(rho*cl_alpha*c*omega^2*R^4);
            Sum(4) = Sum(4) + mb(k);
            dfz(k) = 0.5*rho*c*cl_alpha*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            sz(k) = w(k)*(0.5*(dfz(k)));
            Sum(5) = Sum(5) + sz(k);
            dfx(k) = 0.5*rho*c*cl_alpha*(((Cdo/cl_alpha)*u_t(i)^2  + (u_p(i) * u_t(i)*theta(i)) -u_p(i)^2 ));
            sx(k) = w(k)*0.5*dfx(k);
            Sum(6) = Sum(6) + sx(k);
            dfr(k) = -0.5*rho*c*cl_alpha*beta(i)*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            sr(k) = w(k)*(0.5)*(-beta(i)*dfr(k));
            Sum(7) = Sum(7) + sr(k);
            dfz(k) = 0.5*rho*c*cl_alpha*((u_t(i)^2 * theta(i)) - (u_p(i) * u_t(i)));
            nf(k) = w(k)*0.5*(t(k)*dfz(k));
            Sum(8) = Sum(8) + nf(k);
            dfx(k) = 0.5*rho*c*cl_alpha*(((Cdo/cl_alpha)*u_t(i)^2  + (u_p(i) * u_t(i)*theta(i)) -u_p(i)^2 ));
            nl(k) = w(k)*0.5*t(k)*dfx(k);
            Sum(9) = Sum(9) + nl(k);
        end
        dr_length(j,:) = Sum;
    end
    for a=1:length(Sum);
        blade(i,a) = sum(dr_length(:,a));
    end
    blade(i,10) = 0;
end
r_s_z =blade(:,5);
r_s_x =blade(:,6);
r_s_r =blade(:,7);
m_n_f =blade(:,8);
m_n_l =blade(:,9);
m_n_t =blade(:,10);
%rotating frame
rfx = r_s_x ;
rfy = r_s_r ;
rfz = r_s_z ;
rmx = m_n_f ;
rmy = m_n_t ;
rmz = -m_n_l;

% fixed frame
for i=1:90;
    ffx(i) = (rfy(i)*cos(azm(i)) + rfx(i)*sin(azm(i)))+ (rfy(i+90)*cos(azm(i+90)) + rfx(i+90)*sin(azm(i+90)))+  (rfy(i+180)*cos(azm(i+180)) + rfx(i+180)*sin(azm(i+180))) + (rfy(i+270)*cos(azm(i+270)) + rfx(i+270)*sin(azm(i+270)));
    ffy(i) = (rfy(i)*sin(azm(i)) - rfx(i)*cos(azm(i)))+ (rfy(i+90)*sin(azm(i+90)) - rfx(i+90)*cos(azm(i+90)))+ (rfy(i+180)*sin(azm(i+180)) - rfx(i+180)*cos(azm(i+180))) +(rfy(i+270)*sin(azm(i+270)) - rfx(i+270)*cos(azm(i+270))) ;
    ffz(i) = rfz(i)+rfz(i+90)+rfz(i+180)+rfz(i+270);
    fmx(i) = (rmx(i)*sin(azm(i)) + rmy(i)*cos(azm(i)))+(rmx(i+90)*sin(azm(i+90)) + rmy(i+90)*cos(azm(i+90)))+(rmx(i+180)*sin(azm(i+180)) + rmy(i+180)*cos(azm(i+180))) +(rmx(i+270)*sin(azm(i+270)) + rmy(i+270)*cos(azm(i+270)));
    fmy(i) = (-rmx(i)*cos(azm(i)) + rmy(i)*sin(azm(i)))+(-rmx(i+90)*cos(azm(i+90)) + rmy(i+90)*sin(azm(i+90)))+(-rmx(i+180)*cos(azm(i+180)) + rmy(i+180)*sin(azm(i+180)))+(-rmx(i+270)*cos(azm(i+270)) + rmy(i+270)*sin(azm(i+270)));
    fmz(i) = rmz(i)+rmz(i+90)+rmz(i+180)+rmz(i+270);
end
for i=91:360;
    ffx(i) = ffx(i-90);%(rfy(i)*cos(azm(i)) + rfx(i)*sin(azm(i)))+ (rfy(i+90)*cos(azm(i+90)) + rfx(i+90)*sin(azm(i+90)))+  (rfy(i+180)*cos(azm(i+180)) + rfx(i+180)*sin(azm(i+180))) + (rfy(i+270-360)*cos(azm(i+270-360)) + rfx(i+270-360)*sin(azm(i+270-360)));
    ffy(i) = ffy(i-90);%(rfy(i)*sin(azm(i)) - rfx(i)*cos(azm(i)))+ (rfy(i+90)*sin(azm(i+90)) - rfx(i+90)*cos(azm(i+90)))+ (rfy(i+180)*sin(azm(i+180)) - rfx(i+180)*cos(azm(i+180))) +(rfy(i+270-360)*sin(azm(i+270-360)) - rfx(i+270-360)*cos(azm(i+270-360))) ;
    ffz(i) = ffz(i-90);%rfz(i)+rfx(i+90)+rfx(i+180)+rfx(i+270-360);
    fmx(i) = fmx(i-90);%(rmx(i)*sin(azm(i)) + rmy(i)*cos(azm(i)))+(rmx(i+90)*sin(azm(i+90)) + rmy(i+90)*cos(azm(i+90)))+(rmx(i+180)*sin(azm(i+180)) + rmy(i+180)*cos(azm(i+180))) +(rmx(i+270-360)*sin(azm(i+270-360)) + rmy(i+270-360)*cos(azm(i+270-360)));
    fmy(i) = fmy(i-90);%(-rmx(i)*cos(azm(i)) + rmy(i)*sin(azm(i)))+(-rmx(i+90)*cos(azm(i+90)) + rmy(i+90)*sin(azm(i+90)))+(-rmx(i+180)*cos(azm(i+180)) + rmy(i+180)*sin(azm(i+180)))+(-rmx(i+270-360)*cos(azm(i+270-360)) + rmy(i+270-360)*sin(azm(i+270-360)));
    fmz(i) = fmz(i-90);%rmz(i)+rmz(i+90)+rmz(i+180)+rmz(i+270-360);
end
ndfac_F =1;% rho*pi*R^2*Vtip^2;
ndfac_M =1;% rho*pi*R^3*Vtip^2;

output = [ mean(ffx)/ndfac_F mean(ffy)/ndfac_F mean(ffz)/ndfac_F mean(fmx)/ndfac_M mean(fmy)/ndfac_M mean(fmz)/ndfac_M ];





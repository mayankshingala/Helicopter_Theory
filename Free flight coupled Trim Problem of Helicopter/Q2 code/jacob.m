function J=jacob(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s)
R=Residuefunc(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
new_theta_0=theta_0+0.01*theta_0;
new_R=Residuefunc(new_theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
J1=(new_R-R)/(0.01*new_theta_0);
new_theta_1c=theta_1c+0.01*theta_1c;
new_R=Residuefunc(theta_0, new_theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
J2=(new_R-R)/(0.01*new_theta_1c);
new_theta_1s=theta_1s+0.01*theta_1s;
new_R=Residuefunc(theta_0, theta_1c, new_theta_1s,lambda,mu,phi_s,alpha_s);
J3=(new_R-R)/(0.01*new_theta_1s);
new_lambda=lambda+0.01*lambda;
new_R=Residuefunc(theta_0, theta_1c, theta_1s,new_lambda,mu,phi_s,alpha_s);
J4=(new_R-R)/(0.01*new_lambda);
new_phi_s=phi_s+0.01*phi_s;
new_R=Residuefunc(theta_0, theta_1c, theta_1s,lambda,mu,new_phi_s,alpha_s);
J5=(new_R-R)/(0.01*new_phi_s);
new_alpha_s=alpha_s+0.01*alpha_s;
new_R=Residuefunc(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,new_alpha_s);
J6=(new_R-R)/(0.01*new_alpha_s);
J=[J1 J2 J3 J4 J5 J6];
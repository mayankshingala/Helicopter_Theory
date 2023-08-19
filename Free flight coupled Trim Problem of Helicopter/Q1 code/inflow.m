% function lambda = inflow(lambda_0,mu,alpha_s,CT)
function lambda = inflow(mu,alpha_s,CT)

lambda_h = sqrt(CT/2);         % Initial guess of inflow ratio (typically hover inflow ratio)
% lambda = lambda_h;
lambda_old = lambda_h;          % Old value of Inflow ratio (Initial guess)
maxerrper = 0.0001;             % Maximum error
errper = 1;                     % To initialize the lambda iterative problem 

while errper > maxerrper
    lambda = mu*tan(alpha_s) + CT/(2*sqrt(mu^2 + lambda_old^2));
    errper = abs((lambda - lambda_old)/lambda_old)*100;
    lambda_old = lambda;
end
lambda = lambda_old;
end
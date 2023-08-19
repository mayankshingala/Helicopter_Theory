function lambda=lamdacal(mu,alpha,Ct)
lambdabda_0=mu+Ct/(2*mu);
if mu==0
    lambdabda_0=Ct/2;
end
    error(1)=1;
i=0;
while  error(i+1) > 1e-3    
    f=lambdabda_0-( mu*tan(alpha)+Ct/(2*sqrt(mu^2+lambdabda_0^2)));
    df=1+Ct*lambdabda_0/(2*(lambdabda_0^2+mu^2)^1.5);
    lambda_2=lambdabda_0-f/df;
    i=i+1;
    error(i+1)=abs((lambdabda_0-lambda_2))*100/lambda_2;
    lambdabda_0=lambda_2;
    
end
lambda=lambdabda_0;
  


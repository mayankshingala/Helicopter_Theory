r = r_ct:1:100; 

% No twist: θ_tw= 0 deg
CT_total = sum(CT);
CT_TL_total = sum(CT_TL); 

lambda = sqrt(CT_total/2);
lambda_TL = sqrt(CT_TL_total/2);

theta75_result = ((6*CT_total/(sgm*Cl_a)) + ((3/2)*lambda)) .*ones(1,length(r));
theta75_result_TL = ((6*CT_TL_total/(sgm*Cl_a)) + ((3/2)*lambda_TL)) .*ones(1,length(r));

function [ c, ceq ] = con_evasion( x,param_cost )
%This function summarizes the constraints in the evasion bargaining
%problem
% Recall they must be written as c(x) <= 0
% and if there are NO equality constraints, you must set ceq = []

global tfp lambda 

%Hidden incomes
h = x;

%Total reported income
y_hat = tfp * lambda - f_kappa(h,param_cost)-h;


%Total reported income (=>0) constraint
F1 =  -y_hat;

% %Find the position of time in the vector
% ind_t = find(time_vec == time);
% 
% %Interpolate reported incomes
% e_hat_noscale  = max(10^-5,pchip(tfp*lambda_vec, e_no_evasion(:,ind_t), y_hat));
% w_hat_noscale  = max(10^-5,pchip(tfp*lambda_vec, w_no_evasion(:,ind_t), y_hat));
% 
% %Rescale reported incomes to sum up to total reported income
% e_hat = e_hat_noscale * (y_hat / (e_hat_noscale+w_hat_noscale));
% w_hat = w_hat_noscale * (y_hat / (e_hat_noscale+w_hat_noscale));
% 
% %Find the position of lambda
% ind_lambda = find(lambda_vec == lambda);
% 
% %Employer's should be at least above the threat point 
% % employer's surplus >= 0
% %Worker's should be at least above the threat point
% % workers's surplus >= 0
% if outside_option==1
%     F2 = ce_no_evasion(ind_lambda,ind_t) - (e_hat+he-S_E(w_hat)-T_E(e_hat-S_E(w_hat)));
%     F3 = cw_no_evasion(ind_lambda,ind_t) - (w_hat+hw-S_W(w_hat)-T_W(w_hat));
% else
%     F2 =  - (e_hat+he-S_E(w_hat)-T_E(e_hat-S_E(w_hat)));
%     F3 =  - (w_hat+hw-S_W(w_hat)-T_W(w_hat));
% end
% 
% c_threats  = max(F2,F3);

c = F1;

ceq = [];

end


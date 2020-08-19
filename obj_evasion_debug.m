function [ obj,net_inc,taxes_all,taxes_te,taxes_se,taxes_sw,taxes_tw, e_hat,w_hat,y_hat,e_hat_noscale,w_hat_noscale ] = obj_evasion_debug(x,param_cost)
%This function computes the objective for the Nash bargaining procedure
%when there is tax evasion

global tfp gamma lambda lambda_vec  time time_vec wmin
global e_no_evasion w_no_evasion ce_no_evasion cw_no_evasion
global  display_results wmin_flag

%Hidden incomes
h = x;

%Total reported income
y_hat = tfp * lambda - f_kappa(h,param_cost)-h;
if isreal(y_hat)==0
    warning('y_hat is NOT a real number!')
    keyboard
end

%Find the position of time in the vector
ind_t = find(time_vec == time);

%Find the position of lambda
ind_lambda = find(lambda_vec == lambda);

%Interpolate reported incomes
e_hat_noscale  = max(10^-5,interp1(tfp*lambda_vec, e_no_evasion(:,ind_t), y_hat,'linear','extrap'));
w_hat_noscale  = max(10^-5,interp1(tfp*lambda_vec, w_no_evasion(:,ind_t), y_hat,'linear','extrap'));

%Rescale reported incomes to sum up to total reported income
e_hat = e_hat_noscale * (y_hat / (e_hat_noscale+w_hat_noscale));
w_hat = w_hat_noscale * (y_hat / (e_hat_noscale+w_hat_noscale));
if w_hat < wmin && wmin_flag==1
    w_hat = wmin;
    e_hat = y_hat - w_hat;
end
% Compute total surplus
net_inc = y_hat + h;
taxes_all = max(0,T_E(e_hat-S_E(w_hat)))+max(0,S_E(w_hat))+max(0,S_W(w_hat))+max(0,T_W(w_hat));
taxes_te = max(0,T_E(e_hat-S_E(w_hat)));
taxes_se = max(0,S_E(w_hat));
taxes_sw = max(0,S_W(w_hat));
taxes_tw = max(0,T_W(w_hat));

%surplus_evasion = y_hat + h - T_E(e_hat-S_E(w_hat))-S_E(w_hat)-S_W(w_hat)-T_W(w_hat);
surplus_evasion = y_hat + h - max(0,T_E(e_hat-S_E(w_hat)))-max(0,S_E(w_hat))-max(0,S_W(w_hat))-max(0,T_W(w_hat));
obj = - surplus_evasion;

% If the employer's surplus or the worker's become negative, then obj is
% not a real number! The issue is that fmincon does not always obey the
% nonlinear constraints (it always obeys only the lb<=x<=ub, if algo='sqp' is
% chosen)
% if ~isreal(obj)
    
%pause
% try to fool fmincon :)
if isreal(surplus_evasion)==0  || y_hat<0 
    if display_results==1
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('Warning: objective function in Nash bargaining is NOT a real number!')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
    obj = inf;%[abs(e_surplus_evasion)+abs(w_surplus_evasion)]*100;
    net_inc = NaN;
    taxes_all = NaN;
    taxes_te= NaN;
    taxes_se= NaN;
    taxes_sw= NaN;
    taxes_tw= NaN;
    w_hat   = NaN;
    e_hat   = NaN;
  
    %     else
    %         obj = real(obj)+100000;
end
% end

end
function [ obj ] = obj_no_evasion( e )

% DESCRIPTION:
% This function computes the objective for the Nash bargaining procedure
% when there is no tax evasion

global gamma lambda tfp

w = tfp*lambda - e;

employer_surplus = e-S_E(w)-T_E(e-S_E(w));
worker_surplus   = w -S_W(w)-T_W(w);
obj = - (employer_surplus)^gamma * (worker_surplus)^(1-gamma);

if ~isreal(obj)
    %disp('Warning: objective function in Nash bargaining is NOT a real number')
    obj = abs(real(obj))+1e5; % force fminsearch to stay away from input choices that give bad stuff
end

end


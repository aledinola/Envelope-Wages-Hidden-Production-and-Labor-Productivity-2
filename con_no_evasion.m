function [ c, ceq ] = con_no_evasion( e )

% DESCRIPTION:
% This function summarizes the constraints in the no evasion bargaining problem
% nonlinear inequalities c(x) or equalities ceq(x) defined here. 
% fmincon optimizes such that c(x)<=0 and ceq(x) = 0. 

global lambda tfp 

w = tfp*lambda - e;

c1 = S_E(w) + T_E(e-S_E(w)) - e; % impose employer's surplus >= 0

c2 = S_W(w) + T_W(w) - w; % impose worker's surplus >= 0

c = max(c1, c2);

ceq = [];

end


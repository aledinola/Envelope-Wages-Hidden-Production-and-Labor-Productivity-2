%%% Here we check the results of the Matlab solver using the slow but
%%% robust gridsearch method

L = 300; % grid width
h_vec = linspace(0,(tfp*lambda),L);

% Plot obj func around the "optimal" points found by optimization
% routine
[~,h_ind_solver]=min(abs(x-h_vec));% optimal point according to the nonlinear solver
obj_evasion_matrix = zeros(L,1);
netinc = zeros(L,1);
taxes = zeros(L,1);
taxes_te = zeros(L,1);
taxes_se = zeros(L,1);
taxes_sw = zeros(L,1);
taxes_tw = zeros(L,1);
e_rep = zeros(L,1); % reported profit
w_rep = zeros(L,1); % reported wage
y_rep = zeros(L,1); % reported production

for ie=1:L
    hh = h_vec(ie);
    [obj_evasion_matrix(ie),netinc(ie),taxes(ie),taxes_te(ie),taxes_se(ie),taxes_sw(ie),taxes_tw(ie),e_rep(ie),w_rep(ie),y_rep(ie),e_rep_nos(ie),w_rep_nos(ie)] = obj_evasion_debug(hh,param_cost);
    % if constraints are violated, than obj = INF
end
[~,h_ind] = min(obj_evasion_matrix(:));% optimal point according to the grid search

obj_evasion_matrix(h_ind)
h_evasion_grid(n,t) = h_vec(h_ind);
cost_vec = zeros(length(h_vec),1);
for i=1:length(h_vec)
    cost_vec(i,1) = f_kappa(h_vec(i),param_cost);
end

disp('Optimal points with gridsearch: ')
disp(h_vec(h_ind))
disp('Optimal points with matlab solver: ')
disp(h_vec(h_ind_solver))

figure(1)
subplot(2,1,1)
plot(h_vec,-obj_evasion_matrix,h_vec(h_ind),-obj_evasion_matrix(h_ind),'o',h_vec(h_ind_solver),-obj_evasion_matrix(h_ind_solver),'o','linewidth',2)
legend('obj','grid','solver'), grid on
axis([0 tfp*lambda -inf inf])
xlabel('Hidden Income')
title('Joint Surplus as a function of h')
subplot(2,1,2)
plot(h_vec,netinc,h_vec,taxes,'o','linewidth',2)
legend('y hat + h','taxes'), grid on
axis([0 tfp*lambda -inf inf])

figure(2)
plot(h_vec,cost_vec,h_vec(h_ind),cost_vec(h_ind),'o',h_vec(h_ind_solver),cost_vec(h_ind_solver),'o','linewidth',2)
xlabel('Hidden Income')
legend('cost function','optimal h: grid','optimal h: solver'), grid on
axis([0 tfp*lambda -inf inf])

figure(3)
plot(h_vec,taxes_te,h_vec,taxes_se,h_vec,taxes_sw,h_vec,taxes_tw,'linewidth',2)
xlabel('Hidden Income')
legend('T_E','S_E','S_W','T_W'), grid on
axis([0 tfp*lambda -inf inf])

figure(4)
plot(h_vec,e_rep,h_vec,w_rep,'linewidth',2)
xlabel('Hidden Income')
legend('reported profit','reported wage'), grid on
axis([0 tfp*lambda -inf inf])

figure(5)
plot(h_vec,e_rep_nos,h_vec,w_rep_nos,'linewidth',2)
xlabel('Hidden Income')
legend('reported profit (not rescaled)','reported wage (not rescaled)'), grid on
axis([0 tfp*lambda -inf inf])

figure(6)
plot(e_rep,taxes_te,'linewidth',2)
xlabel('Reported Profit')
legend('T_E'), grid on
axis([0 tfp*lambda -inf inf])

figure(7)
plot(w_rep,taxes_se,w_rep,taxes_sw,w_rep,taxes_tw,'linewidth',2)
xlabel('Reported Wage')
legend('S_E','S_W','T_W'), grid on
axis([0 tfp*lambda -inf inf])

figure(8)
plot(h_vec,y_rep,'linewidth',2)
xlabel('Hidden Income')
legend('y hat: reported production'), grid on
axis([0 tfp*lambda -inf inf])


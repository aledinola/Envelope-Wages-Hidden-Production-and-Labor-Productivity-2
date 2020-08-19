% check some results
tt = 1; % pick year (from 1 to 15)
nn_lambda = 7; %n_lambda; % from 1 to 10
%keyboard
% check wage
figure(1)
plot(1:nn_lambda,w_no_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
legend('w NO evasion','min wage')
figure(2)
plot(1:nn_lambda,w_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
legend('w evasion','min wage')
figure(3)
plot(tfp*lambda_vec(1:nn_lambda),w_no_evasion(1:nn_lambda,tt),'-o',y_evasion(1:nn_lambda,tt),w_hat_noscale(1:nn_lambda,tt),'linewidth',2),grid on
legend('w NO evasion','w evasion (interp)')
% check profit
figure(4)
plot(1:nn_lambda,e_no_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
legend('e NO evasion','min wage')
figure(5)
plot(1:nn_lambda,e_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
legend('e evasion','min wage')
figure(6)
plot(tfp*lambda_vec(1:nn_lambda),e_no_evasion(1:nn_lambda,tt),'-o',y_evasion(1:nn_lambda,tt),e_hat_noscale(1:nn_lambda,tt),'linewidth',2),grid on
legend('e NO evasion','e evasion (interp)')

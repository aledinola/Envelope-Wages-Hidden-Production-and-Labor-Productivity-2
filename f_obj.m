function [ objfun_smm ] = f_obj( param_vec )

%% Description of the function
% INPUTS:
%  param_vec: parameters to be estimated
%  It contains 15 z(t)'s values
% OUTPUT: 
% objfun_smm = Squared relative distance b/w model targets and data targets

%% DECLARE GLOBALS:

global gamma time lambda tfp wmin time_vec LB UB
global cost_specif display_results
global lambda_vec w_no_evasion e_no_evasion cw_no_evasion ce_no_evasion
global cubic interpol outside_option
global S_E_year S_W_year T_W_year year_vec  %#ok<NUSED>
global experiment experiment_year do_estimation 
global  labor_prod_data av_wages informal_size wmin_vec start_year  final_year reform_year
global n t h_evasion est_algo wmin_flag la freeze_productivity do_experiments experiment_num

%% SET SOME CONTROL FLAGS:

% Set here the algorithm for constrained minimization:
%---------------------- NO TAX EVASION -----------------------------------% 
% N.B. it seems it is better to use fmincon
algo_noev = 'fmincon_algo'; % fmincon is a derivative-based algorithm, fast but maybe doesn't
                       % handle well non-smooth functions, non convex
                       % constraints etc.
%algo_noev = 'fminsearchcon_algo'; % downloaded from mathworks. It is a modified fminsearch routine 
%                                   % that can handle nonlinear constraints
                        
%--------------------- TAX EVASION ECONOMY -------------------------------% 
%(better use fminsearchcon)                    
% algo = 'fmincon_algo'; % fmincon is a derivative-based algorithm, fast but maybe doesn't
%                        % handle well non-smooth functions, non convex
%                        % constraints etc.
% algo = 'fminsearch_algo'; % fminsearch with constraints directly imposed in the obj function
%                          % infinite penalization method: if any constraint
%                          % is violated, set obj = Inf. The pros of
%                          % fminsearch is that it is a derivative-free
%                          % algorithm
algo = 'fminsearchcon_algo'; % downloaded from mathworks. It is a modified fminsearch routine 
                        % that can handle nonlinear constraints
% algo = 'genetic_algo'; % most robust method, global 
%-------------------------------------------------------------------------%

    
display_results = 0; % = 1 if u wanna display flags when running the code (for debug mode only)
grid_search     = 0; % set = 1 if you want to check results using grid search
%plot_results    = 0; % = 1 if you wanna plot results and compute statistics
    
interpol = 1; % choose whether to interpolate the taxes or to use the estimated log specification
              % (the latter uses the old code)
outside_option = 0; % set = 0 if the outside option in the tax evasion economy is zero
                    % set = 1 if the outside option in the tax evasion
                    % economy is the after-tax income if NO evasion
n_grid_init = 100; % % width of grid to choose initial conditions (refinement)              
experiment = 0; % set = 1 if u wanna perform counterfactual experiments
experiment_year = 2000; 

%% SET PARAMETER VALUES (TO BE ESTIMATED)
% check bounds (if in estimation mode)

if do_estimation == 1
    param_vec_col = size(param_vec,1);
    LB_col = size(LB,1);
    if param_vec_col~=LB_col
        param_vec = param_vec';
    end
    if isempty(find(param_vec<LB))==0 %#ok<EFIND>
        error('Bounds are violated!')
    elseif isempty(find(param_vec>UB))==0 %#ok<EFIND>
        error('Bounds are violated!')
    end
end

if strcmp(cost_specif,'c6') % DOUBLE THRESHOLD
    %guess = [alpha; beta; theta; gamma; var_lambda; scale; delta; x_lowbar; x_highbar];
    alpha = param_vec(1);
    beta = param_vec(2);
    theta = param_vec(3);
    gamma = param_vec(4);
    var_lambda = param_vec(5);
    scale = param_vec(6); %#ok<*NASGU>
    delta = param_vec(7);
    x_lowbar = param_vec(8);
    x_highbar = param_vec(9);
    %PASS PARAMETERS FOR COST FUNCTION
    param_cost.alpha = alpha;
    param_cost.beta = beta;
    param_cost.theta = theta;
    
elseif strcmp(cost_specif,'c5') % LOGISTIC
    %guess = [alpha; beta; theta; xbar; gamma; var_lambda; scale; delta];
    alpha = param_vec(1);
    beta = param_vec(2);
    theta = param_vec(3);
    gamma = param_vec(4);
    var_lambda = param_vec(5);
    scale = param_vec(6);
    delta = param_vec(7);
    %PASS PARAMETERS FOR COST FUNCTION
    param_cost.alpha = alpha;
    param_cost.beta = beta;
    param_cost.theta = theta;

elseif strcmp(cost_specif,'c7') % EXPONENTIAL
   
    beta = param_vec(1);
    theta = param_vec(2);
    gamma = param_vec(3);
    var_lambda = param_vec(4);
    tfp_vec = param_vec(5:end);
    %PASS PARAMETERS FOR COST FUNCTION
    param_cost.beta = beta;
    param_cost.theta = theta;
    
else
    %guess = [alpha; beta; theta; gamma; var_lambda; scale; delta];
    alpha = param_vec(1);
    beta = param_vec(2);
    theta = param_vec(3);
    gamma = param_vec(4);
    var_lambda = param_vec(5);
    scale = param_vec(6);
    delta = param_vec(7);
    %PASS PARAMETERS FOR COST FUNCTION
    param_cost.alpha = alpha;
    param_cost.beta = beta;
    param_cost.theta = theta;
end

%% SET THE GRID FOR STATE VARIABLES

n_lambda = 10;     %Number of grid points for ability
[loglambda_vec,pi_initial] = tauchen(n_lambda,-var_lambda/2,0,var_lambda^0.5,5.0);
%[loglambda_vec,pi_initial]=tauchenhussey(n_lambda,0,0,var_lambda^0.5,var_lambda^0.5);
lambda_prob = pi_initial(1,:)';
lambda_vec1 = exp(loglambda_vec);
lambda_vec = lambda_vec1/(lambda_vec1'*lambda_prob);
% check distribution of lambdas
% plot(lambda_vec,lambda_prob)

%% SET TIME AND LOAD THE TFP SERIES

t_time = final_year-start_year+1; %number of time periods
time_vec = transpose(start_year:final_year);
%tfp_vec = tfp_vec1*scale;


if do_experiments==1 && freeze_productivity==1
    % set TFP always equal to 2000 level
    tfp_vec = ones(t_time,1)*tfp_vec(1);
end


try
if abs(length(tfp_vec) - t_time) > 0
    error('The TFP series do not correspond to the time period')
end
catch ME
    keyboard
end

% LOAD THE TAXES SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Tax_functions.mat')

%%%%%%%%% INTERPOLATION?
if interpol==1
    % CHOOSE PREFERRED INTERP METHOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cubic = 0; % if cubic=0, then do linear interp
               % it seems cubic creates problems with fmincon
end


%% COMPUTE THE NO-EVASION ECONOMY AND THE THREAT POINTS

%Track time
tic

%Print
disp('=================================')
disp('Number of gridpoints for lambda:')
disp(n_lambda)
disp(' ')
disp('Computing the no-evasion problem.')
disp('Please wait...')

% Initialize the storing matrices
e_no_evasion  = zeros(n_lambda, t_time); % employers income
w_no_evasion  = zeros(n_lambda, t_time); % workers income
ce_no_evasion = zeros(n_lambda, t_time); % employers after-tax income
cw_no_evasion = zeros(n_lambda, t_time); % workers aftre-tax income

% Control matrices
fval_no_evasion     = zeros(n_lambda, t_time); % objective function at optimum
exitflag_no_evasion = zeros(n_lambda, t_time); % exiflag of the minimization
inactive  = zeros(n_lambda, t_time); % for each t,n is the plant active or not?
                                     % 1 = inactive, 0 = active

%Options
switch algo_noev
    case 'fmincon_algo'
        options=optimset('MaxFunEvals',3000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off','Algorithm','sqp');
    case 'fminsearchcon_algo'
        options=optimset('MaxFunEvals',3000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off');
end


for t = 1:t_time
    
    %Set the time and the tfp
    time = time_vec(t);
    tfp = tfp_vec(t);
    wmin = wmin_vec(t);
    guess = gamma*tfp*lambda_vec(1); % analytical solution of NB for e = profit
    fprintf('Year is %d\n',time)
    
    for n = 1:n_lambda
        %Set the ability
        lambda = lambda_vec(n);
        % If TFP*lambda < minimum wage, then absent production pair cannot profitably operate
        if tfp*lambda>=wmin
            %Solve the problem
            switch algo_noev
                case 'fmincon_algo'
                    [x,fval,exitflag] = fmincon('obj_no_evasion',guess,...
                        [],[],[],[],0,tfp*lambda,'con_no_evasion',options);
                    % x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options) subjects the minimization to the nonlinear
                    % inequalities c(x) or equalities ceq(x) defined in nonlcon. fmincon optimizes such that c(x) <= 0 and ceq(x) = 0.
                    % If no bounds exist, set lb = [] and/or ub = []
                case 'fminsearchcon_algo'
                    [x,fval,exitflag] = fminsearchcon('obj_no_evasion',guess,0,tfp*lambda,[],[],'con_no_evasion',options);
                    % usage: x=FMINSEARCHCON(fun,x0,LB,UB,A,b,nonlcon,options)
            end
            
            %Write the results
            
            e_no_evasion(n,t) = x;
            w_no_evasion(n,t) = tfp*lambda - x;
            
            fval_no_evasion(n,t) = -fval;
            guess = x;

            %Exit if flag is not OK
            exitflag_no_evasion(n,t) = exitflag;

            if exitflag < 1
                warning('The NO-evasion problem has no feasible solution for gridpoint %d\n',n)
            end

            %Check whether the constraint is satisfied
            if con_no_evasion( x ) > 10^-10
                if do_estimation==1
                    objfun_smm = 10^100000;
                    return
                else
                    %keyboard
                    error('The no-evasion problem is not solved properly. Check constraints.')
                end
            end

            %Check whether e and w are above or eq to zero
            if e_no_evasion(n,t) < -10^-10 || w_no_evasion(n,t) < -10^-10
                if do_estimation==1
                    objfun_smm = 10^100000;
                    return
                else
                    error('Incomes are negative. Check the solution.')
                end
            end
            
            % Impose the min.wage constraint
            if w_no_evasion(n,t)<wmin && wmin_flag == 1
                disp('MIN WAGE IS BINDING')
                w_no_evasion(n,t) = wmin;
                e_no_evasion(n,t) = tfp*lambda - wmin;
            end

            %After-tax incomes
            ce_no_evasion(n,t) = e_no_evasion(n,t) - S_E(w_no_evasion(n,t)) - T_E(e_no_evasion(n,t)-S_E(w_no_evasion(n,t)));
            cw_no_evasion(n,t) = w_no_evasion(n,t) - S_W(w_no_evasion(n,t)) - T_W(w_no_evasion(n,t));
        
        else % case when y < min.wage ==> inactive flag is set to 1 
            fprintf('Firm with lambda %d is INACTIVE\n', n )
            e_no_evasion(n,t)  = 0;
            w_no_evasion(n,t)  = 0;
            ce_no_evasion(n,t) = 0;
            cw_no_evasion(n,t) = 0;
            inactive(n,t)      = 1;
        end % end check for inactivity

    end % end loop over n
end % end loop over t

la_vec = nan(t_time,1);
for t=1:t_time
% Collect here the first active lambda
    la_vec(t) = find(inactive(:,t)==0, 1 );
end

%Print
disp('The no-evasion problem is computed!')

%==========================================================================%
% %%% Some plots to check
% tt = 1; % pick year (from 1 to 15)
% nn_lambda = 7; %n_lambda; % from 1 to 10
% % check wage
% figure
% plot(1:nn_lambda,w_no_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
% legend('w NO evasion','min wage'),title('WAGE')
% % check profit
% figure
% plot(1:nn_lambda,e_no_evasion(1:nn_lambda,tt),'-o',1:nn_lambda,wmin_vec(tt)*ones(nn_lambda,1),'linewidth',2),grid on
% legend('e NO evasion','min wage'),title('PROFIT')
%==========================================================================%


%% HOW TO DEAL WITH INACTIVITY

% Suppose for a given year t all pairs s.t. lambda<=lambda_t are inactive
% This happens if z(t)*lambda(t) < wmin(t)
% When solving the economy with tax evasion, we should consider only the
% lambdas s.t. lambda > lambda_t
% Then we have to renormalize the probability density when computing the
% aggregate variables

%% COMPUTE THE TAX-EVASION ECONOMY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Track time
% tic

%Print
disp('=================================')
disp('Computing the tax evasion problem.')
disp('Please wait...')

% Initiate the storing matrices
% Wrt the inefficient model, we have only ONE choice variable for tax
% evasion
h_evasion = zeros(n_lambda, t_time); % hidden income
if grid_search==1
    h_evasion_grid = zeros(n_lambda, t_time); % to check results
end
h_evasion_prod = zeros(n_lambda, t_time);
y_prod         = zeros(n_lambda, t_time);
he_evasion     = zeros(n_lambda, t_time);
hw_evasion     = zeros(n_lambda, t_time);
e_hat_noscale  = zeros(n_lambda, t_time); %employers reported income
w_hat_noscale  = zeros(n_lambda, t_time); %workers reported income
e_evasion      = zeros(n_lambda, t_time); %employers reported income
w_evasion      = zeros(n_lambda, t_time); %workers reported income
y_evasion      = zeros(n_lambda, t_time); %total reported income
cw_evasion      = zeros(n_lambda, t_time); %workers reported income (AFTER TAX)
%e_surplus_evasion = zeros(n_lambda, t_time);
%w_surplus_evasion = zeros(n_lambda, t_time);

% Control matrices
fval_evasion = zeros(n_lambda, t_time); % objective function at optimum
exitflag_evasion = zeros(n_lambda, t_time); % exiflag of the minimization

% Options
switch algo
    case 'fmincon_algo'
        options=optimset('MaxFunEvals',4000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off','Algorithm','sqp');
    case 'fminsearch_algo'
        options=optimset('MaxFunEvals',10000,'MaxIter',10000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off');
    case 'fminsearchcon_algo'
        options=optimset('MaxFunEvals',10000,'MaxIter',10000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off');
    case 'genetic_algo'
%         options=gaoptimset('MaxFunEvals',4000,'FunValCheck','off','TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10,'disp','off','Algorithm','sqp');
end
    
for t = 1:t_time
    
    %Set the time and the tfp
    time = time_vec(t);
    fprintf('Year is %d\n',time)
    tfp = tfp_vec(t);
    wmin = wmin_vec(t);
    lambda = lambda_vec(1);
    la = la_vec(t); % First lambda where firm is ACTIVE
     %gridsearch % for each year t, find a good initial condition
    
     for n = la:n_lambda
        fprintf('Gridpoint no. %d\n',n) 

    %Set the ability
    lambda = lambda_vec(n); 
    % check if initial condition satisfies constraints; if not, do again
    gridsearch %to find a suitable guess
    if con_evasion(guess,param_cost) > 0
        gridsearch
    end
    
    %Solve the problem
        switch algo
            case 'fmincon_algo'
                [x,fval,exitflag] = fmincon(@(xx) obj_evasion(xx,param_cost),guess,...
                           [],[],[],[],[0],[tfp*lambda],@(x) con_evasion(x,param_cost),options);
            case 'fminsearch_algo'
                 [x,fval,exitflag] = fminsearch(@(xx) obj_evasion(xx,param_cost),guess,options);
            case 'fminsearchcon_algo'
                 [x,fval,exitflag]=fminsearchcon(@(xx) obj_evasion(xx,param_cost),guess,[0],[tfp*lambda],[],[],@(x) con_evasion(x,param_cost),options);
                 % usage: x=FMINSEARCHCON(fun,x0,LB,UB,A,b,nonlcon,options)
            case 'genetic_algo'
                 [x,fval,exitflag] = ga(@(xx) obj_evasion(xx,param_cost),2,[],[],[],[],[0],[tfp*lambda],@con_evasion);
        end
        
        % Grid search to check results (for given year t and lambda n)
        if grid_search==1
            taxev_debug
        end
        %%%
                       
        fval_evasion(n,t) = -fval;
        exitflag_evasion(n,t) = exitflag;
        %guess = x;        
        %Exit if flag is not OK
        if exitflag < 1
            warning('The evasion problem has no feasible solution for gridpoint %d\n',n)
        end
        
        %Check whether the constraints are satisfied
        if con_evasion(x,param_cost) > 10^-10
            error('The evasion problem is not solved properly. Check constraints.')
        end               

        %Write the results
        h_evasion(n,t) = x;
        y_prod(n,t) = (tfp_vec(t)*lambda_vec(n));
        h_evasion_prod(n,t) = h_evasion(n,t)./y_prod(n,t);
        he_evasion(n,t) = gamma*x; % non-reported (hidden) profit
        hw_evasion(n,t) = (1-gamma)*x; % non-reported wage
        
        %Total reported income
        y_evasion(n,t) = tfp * lambda - f_kappa(he_evasion(n,t)+hw_evasion(n,t),param_cost) ...
                         - he_evasion(n,t) - hw_evasion(n,t);
%         y_evasion(n,t)=tfp*lambda-he_evasion(n,t)-hw_evasion(n,t);
        %Interpolate reported incomes
%         e_hat_noscale  = max(10^-5,pchip(tfp*lambda_vec, e_no_evasion(:,t), y_evasion(n,t)));
%         w_hat_noscale  = max(10^-5,pchip(tfp*lambda_vec, w_no_evasion(:,t), y_evasion(n,t)));
        e_hat_noscale(n,t)  = max(10^-5,interp1(tfp*lambda_vec(la:end), e_no_evasion(la:end,t), y_evasion(n,t),'linear','extrap'));
        w_hat_noscale(n,t)  = max(10^-5,interp1(tfp*lambda_vec(la:end), w_no_evasion(la:end,t), y_evasion(n,t),'linear','extrap'));
        
        %Rescale reported incomes to sum up to total reported income
        e_evasion(n,t) = e_hat_noscale(n,t) * (y_evasion(n,t) / (e_hat_noscale(n,t)+w_hat_noscale(n,t)));
        w_evasion(n,t) = w_hat_noscale(n,t) * (y_evasion(n,t) / (e_hat_noscale(n,t)+w_hat_noscale(n,t)));
        
        % After-tax reported incomes (profit and wages)
        cw_evasion(n,t) = w_evasion(n,t) - S_W(w_evasion(n,t)) - T_W(w_evasion(n,t));
                       
    end  % end loop over n
    
end % end loop over t

%Print
disp('The evasion problem is computed!')

%Track time
toc

% Display figures only if we are in no-estimation mode
if do_estimation==0
    figures_debug
    close all
end


%% COMPUTE TARGETS
% clear all; clc
% load('temp_results.mat')
% (1) Size of informal economy t=1,..,15
% (2) Average (observed) wage t=1,..,15
% (3) GDP per employee t=1,..,15
% (4) Gini coefficient 
% Compute simulated moments 



GDP_year     = tfp_vec;         % GDP per employee: trivial, in aggregate Y = z (tfp)
H_year       = zeros(t_time,1); % Aggregate hidden income he+hw, for each year
H_year_check = zeros(t_time,1);
w_evasion_year = zeros(t_time,1); % average reported wages
e_evasion_year = zeros(t_time,1); % average reported profit
y_evasion_year = zeros(t_time,1); % reported production Y hat (two ways to compute it)
Gini_w_evasion_year = zeros(t_time,1); % Gini on reported wage, by year 
Gini_cw_evasion_year = zeros(t_time,1); % Gini on reported wage, after tax, by year 
Gini_y_evasion_year = zeros(t_time,1); % Gini on reported total income, by year 
L_year = zeros(t_time,1); % aggregate efficiency loss
T_E_aggregate = zeros(t_time,1); % aggregate corporate taxes
T_W_aggregate = zeros(t_time,1); % aggregate labor income taxes
T_year        = zeros(t_time,1); % Tax revenues raised by gov.
wmin_year = zeros(t_time,1); % share of production pairs at the minimum wage
inactive_year = zeros(t_time,1); % share of inactive production pairs
av_lambda= sum(lambda_prob.*lambda_vec);
%lambda_prob=lambda_prob/av_lambda;

for t= 1:t_time
    for n=1:n_lambda
        inactive_year(t) = inactive_year(t) + inactive(n,t)*lambda_prob(n);
    end % inactive is a dummy equal to 1 if INACTIVE, 0 if ACTIVE
end

for t= 1:t_time
    time = time_vec(t);
    tfp = tfp_vec(t);
    la = la_vec(t); % n = la is the lowest active pair
    lambda_prob_active = lambda_prob/(1-inactive_year(t)); % vector (n_lambda-la+1)*1 or n_lambda*1
    disp(['Sum of prob mass for active firms: ',num2str(sum(lambda_prob_active(la:end)))])
    fprintf('E(lambda|lambda>=lambda_min): %f\n', sum(lambda_vec(la:end).*lambda_prob_active(la:end)))
    % Gini on reported labor income
    [Gini_w_evasion_year(t),~]=gini(lambda_prob(la:end),max(w_evasion(la:end,t),0));
    % Gini on total reported income
    [Gini_y_evasion_year(t),~]=gini(lambda_prob(la:end),max(y_evasion(la:end,t),0));
    % Gini coefficient for after-tax labor income (disposable income) 
    [Gini_cw_evasion_year(t),~]=gini(lambda_prob(la:end),max(cw_evasion(la:end,t),0));
    
    for n=la:n_lambda % loop only over active production pairs
        
        % TO DO: I should normalize the probability(lambda) by the share of
        % active firms
        lambda = lambda_vec(n);
        H_year(t) = H_year(t)+ he_evasion(n,t)*lambda_prob_active(n)+ hw_evasion(n,t)*lambda_prob_active(n);
        H_year_check(t) = H_year_check(t)+ h_evasion(n,t).*lambda_prob_active(n);
        L_year(t) = L_year(t)+ f_kappa(h_evasion(n,t),param_cost)*lambda_prob_active(n);
        % f_kapp can return a vector if input is vector
        w_evasion_year(t) = w_evasion_year(t)+ w_evasion(n,t).*lambda_prob_active(n);
        e_evasion_year(t) = e_evasion_year(t)+ e_evasion(n,t).*lambda_prob_active(n);
        y_evasion_year(t) = y_evasion_year(t)+ y_evasion(n,t).*lambda_prob_active(n);
        d = w_evasion(n,t)<=wmin_vec(t)*1.0001;
        wmin_year(t)      = wmin_year(t) + d * lambda_prob_active(n);
        
        T_E_aggregate(t) = 0;
        T_W_aggregate(t) = 0;
        T_E_aggregate(t) = T_E_aggregate(t) + T_E(e_evasion(n,t))*lambda_prob_active(n);
        T_W_aggregate(t) = T_W_aggregate(t) + T_W(w_evasion(n,t))*lambda_prob_active(n);
        T_year(t) = T_W_aggregate(t) + T_E_aggregate(t);
    end
end

% reported production, by year
Y_hat_year = GDP_year - H_year - L_year; % resource constraint of the economy
informal_size_model = H_year./GDP_year;
informal_size_model_reported = H_year./Y_hat_year;

reform_year_index = reform_year-2000+1;
if t_time==1
    Gini_w_evasion_pre = Gini_w_evasion_year;
    Gini_w_evasion_post = Gini_w_evasion_year;
    Gini_y_evasion_pre = Gini_y_evasion_year;
    Gini_y_evasion_post = Gini_y_evasion_year;
else  % Gini on labor income and Gini on total income
    Gini_w_evasion_pre = mean(Gini_w_evasion_year(1:reform_year_index));
    Gini_w_evasion_post = mean(Gini_w_evasion_year(reform_year_index+1:end));
    Gini_y_evasion_pre = mean(Gini_y_evasion_year(1:reform_year_index));
    Gini_y_evasion_post = mean(Gini_y_evasion_year(reform_year_index+1:end));
    Gini_y_evasion = mean(Gini_y_evasion_year(1:end));
    Gini_cw_evasion = mean(Gini_cw_evasion_year(1:end));
    Gini_cw_evasion_pre = mean(Gini_cw_evasion_year(1:reform_year_index));
    Gini_cw_evasion_post = mean(Gini_cw_evasion_year(reform_year_index+1:end));
end

% Additional Target: Minimum wage share
wmin_share_data = 0.33;
wmin_share_model = mean(wmin_year);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMM (objfun_smm will become the output of a function called 'model_smm')

if t_time>1
model_moments = [informal_size_model_reported;
                 w_evasion_year; % average reported wages
                 Y_hat_year; % reported total production
                 Gini_cw_evasion;%Gini_y_evasion;
                 wmin_share_model];
% gini_pre_data = 0.312160966259722;
% gini_post_data = 0.346931640261111;
gini_data = 0.332;

data_moments = [informal_size;
                av_wages; % average wages obs in the data
                labor_prod_data; % GDP per employee, observed in data
                gini_data;
                wmin_share_data];

%data_moments = informal_size;
weights = ones(length(model_moments),1);
weights(1:15) = 2; % importance given to informal size
distance = weights.*(abs(data_moments-model_moments)./data_moments);
objfun_smm = sum(distance.^2);

% save estimation_temp_results.mat
% save paramfile alpha beta theta gamma var_lambda

%Display intermediate results on the screen

fprintf('=========================================\n')
fprintf('MODEL FIT SO FAR\n')
disp('Informal sector size, data vs model')
disp([informal_size,informal_size_model_reported])
fprintf('-----------------------------------------\n')
disp('Reported Wages, data vs model')
disp([av_wages,w_evasion_year])
fprintf('------------------------------------------\n')
disp('Reported Labor prod, data vs model')
disp([labor_prod_data,Y_hat_year])
fprintf('-------------------------------------------\n')


fprintf('=============================\n')
fprintf('PARAMETERS, current values\n')

fprintf('beta(level) = %.12f\n',beta)
fprintf('Theta(convexity) = %.12f\n',theta)
fprintf('Gamma (Employer barg. power) = %.6f\n',gamma)
fprintf('Var.lambda shocks = %.6f\n',var_lambda)
fprintf('tfp %.6f\n',tfp_vec)
fprintf('=============================\n')
fprintf('Obj fun = %.4f\n',objfun_smm)
fprintf('=============================\n')


% save intermediate results in a txt file
if do_estimation==1
    fid = fopen('paramfile_results.txt','at');
    fprintf(fid,'=============================\n');
    fprintf(fid,'minim. routine: %s\n',est_algo);
    fprintf(fid,'PARAMETERS, current values\n');
    fprintf(fid,'beta = %.12f\n',beta);
    fprintf(fid,'theta = %.12f\n',theta);
    fprintf(fid,'gamma = %.6f\n',gamma);
    fprintf(fid,'var_lambda = %.6f\n',var_lambda);
    fprintf(fid,'tfp = %.6f\n',tfp_vec);
    fprintf(fid,'=============================\n');
    fprintf(fid,'Obj fun = %.4f\n',objfun_smm);
    fprintf(fid,'=============================\n');
    fclose(fid);
end    
    
end

% Display figures only if we are in no-estimation mode

if do_estimation==0
    
    if do_experiments==0
        make_plots
        %make_plots_sensitivity
    elseif do_experiments==1
        % change here manually!
        if experiment_num == 1
            informal_z_fixed = informal_size_model_reported;
            save informal_z_fixed informal_z_fixed;
        elseif experiment_num==2
            informal_alltax_fixed = informal_size_model_reported;
            save informal_alltax_fixed informal_alltax_fixed;
        elseif experiment_num==3
            informal_tw_only = informal_size_model_reported;
            save informal_tw_only informal_tw_only;
        elseif experiment_num==4
            informal_te_only = informal_size_model_reported;
            save informal_te_only informal_te_only
        elseif experiment_num==5
            informal_s_only = informal_size_model_reported;
            save informal_s_only informal_s_only
            
            
        end
    end
    
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   "Envelope Wages, Hidden Production and Labor Productivity"            %
%                                                                         %
%       by A. Di Nola, G. Kocharkov and A. Vasilev                        % 
%                                                                         %
%       published on the B.E. Journal of Macroeconomics (2019)            %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
clear; clc; close all
format long g

disp('Envelope Wages, hidden production and labor productivity')
disp('by A.Di Nola, G.Kocharkov and A.Vasilev')
disp('published on the B.E. Journal of Macroeconomics 2019 (Advances)')
fprintf('\n')

global cost_specif do_estimation year_vec labor_prod_data av_wages informal_size 
global wmin_vec start_year  final_year reform_year
global LB UB est_algo wmin_flag 
global freeze_te freeze_se freeze_sw freeze_tw freeze_productivity experiment_year do_experiments experiment_num


%% Set some flags
% If you want to simply run the model at the calibrated values, 
% in order to replicate the resuts in the paper,
% please set do_estimation=0, do_experiments=0 and param_from_file=1.

param_from_file = 1;     % Choose parameters from .mat file
do_estimation   = 0;     % if 0, then simply run the model at x = guess  
do_parallel     = 0;     % set it = 1 if parallel computing
wmin_flag       = 1;     % set it = 1 if wanna impose the min.wage constraint in the 
                         % NO-tax-evasion problem
                         
%% Set some extra flags for the counterfactual experiments

%-------------- 5 experiments----------------------------------------------

% 1 - keep productivity fixed at 2000 level
% 2 - keep all taxes fixed
% 3 - change only T_W (other taxes and prod fixed)
% 4 - change only T_E 
% 5 - change only S_E & S_W 
%--------------------------------------------------------------------------

do_experiments = 0; % 1 = do experiments
experiment_num = 5; % we have 5 experiments so far
% If, for ex., freeze_te = 1, then the model is solved keeping the
% corporate business tax "freezed" at its 2000 level.
freeze_te = 1; % corporate business tax
freeze_se = 0; % employers social contribution
freeze_sw = 0; % workers social contribution
freeze_tw = 1; % personal income tax

freeze_all_taxes = 0; % set = 1 if u wanna perform counterfactual experiments
if freeze_all_taxes==1
   freeze_te = 1; % corporate business tax
   freeze_se = 1; % employers social contribution
   freeze_sw = 1; % workers social contribution
   freeze_tw = 1; % personal income tax
end
experiment_year = 2000; 
freeze_productivity = 1;


%% Flags for the minimization routine

% Choose estimation algo for the minimization routine (5 options)
est_algo = 'nlopt';
% est_algo = 'simulan';
% est_algo = 'ga';
% est_algo = 'simulan_mathworks';
% est_algo = 'fminsearch';

%% Flags for cost function specification

% cost_specif = 'c1'; % alpha+beta*h^theta + delta(h/y-h_past/y_past)^2
% cost_specif = 'c2'; % % alpha+beta*h^theta + delta(h/y-h_past/y_past)^2
% cost_specif = 'c3'; % cost depends on h/y ==> high prod evade more
% cost_specif = 'c4'; % cost  = beta * (alpha*h)^theta/(tfp_vec(t)^theta);
% cost_specif = 'c5'; % LOGISTIC
% cost_specif = 'c6'; % DOUBLE THRESHOLD
cost_specif = 'c7'; % EXPONENTIAL (best specification so far)

%% Load fixed/external parameters

% Fixed parameters (they do not change in the estimation)
% Time parameters
start_year = 2000;
final_year = 2014;
reform_year = 2008;
if start_year<2000 || final_year>2014
    error('Time range is not feasible!')
end
load data
year_vec = data(:,1);
labor_prod_data = data((start_year:final_year)-2000+1,2); % this is real GDP per employee
% based on observed GDP
tfp_vec1 = data((start_year:final_year)-2000+1,2); % these are z(t)
av_wages = data(:,3); % average reported wages
% informal sector size (for 2014 is missing!!)
informal_size = data(:,5)/100; % in excel file it is in %
wmin_vec = data(:,4); % minimum wages

%% SET PARAMETER VALUES:
% These values are used as initial points for the minimization if
% do_estimation = 1. Otherwise we run f_obj(x) where x = parameters

if strcmp(cost_specif,'c1')
    % c1 cost specification
    alpha = 0;%0.927020;
    beta = 0.011;%0.011507;
    theta = 1.4;%33008;
    gamma = 0.95;%0.080461;
    var_lambda = 1.127645;
    scale  = 1.4;%0.101393;
    %delta = 100000;
    delta = 0;
elseif strcmp(cost_specif,'c2')
    % quadratic but we don't match labor prod
    alpha = 0.000000;
    beta = 0.010060;
    theta = 1.533270;
    delta = 0.000000;
    gamma  = 0.257572;
    var_lambda  = 0.931397;
    scale  = 0.154671;
    
elseif strcmp(cost_specif,'c4')
    alpha = 6.5;%0.927020;
    beta = 500;%0.011507;
    theta = 1.4;%33008;
    gamma = 0.95;%0.080461;
    var_lambda = 0.15;%1.127645;
    scale  = 1.4;%0.101393;
    %delta = 100000;
    delta = 0;%100000;
    
elseif strcmp(cost_specif,'c5') % LOGISTIC
    % alpha: related to the asymptotical value (max as h==>inf)
    % beta:  inflexion point
    % theta: steepness of the logistic
    alpha = 7507.984257;
    beta = 4100; %3988.150639; both are good!
    theta = 0.024222;
    gamma = 0.95;%0.080461;
    var_lambda = 0.15;%1.127645;
    scale  = 1.4;%0.101393;
    delta = 0;%100000;
    
elseif strcmp(cost_specif,'c6') % DOUBLE THRESHOLD
    alpha = 1;%0.927020;
    beta = 0.0002;
    theta = 2;
    x_lowbar = 0.1;
    x_highbar = 0.51; % 0.5
    gamma = 0.95;
    var_lambda = 0.15;
    scale  = 1.4;
    delta = 0;
    
elseif strcmp(cost_specif,'c7') % EXPONENTIAL
    
    beta = 0.000000022441;
    theta = 0.006971030922;
    gamma = 0.975955;
    var_lambda = 0.486333;
    tfp_guess = zeros(15,1);
    tfp_guess(1) = 11314.998723;
    tfp_guess(2) = 11103.999171;
    tfp_guess(3) = 10978.798849;
    tfp_guess(4) = 11143.214828;
    tfp_guess(5) = 11194.796602;
    tfp_guess(6) = 11266.737811;
    tfp_guess(7) = 11302.129582;
    tfp_guess(8) = 11368.968917;
    tfp_guess(9) = 11696.609665;
    tfp_guess(10) = 11512.356740;
    tfp_guess(11) = 11737.187976;
    tfp_guess(12) = 11971.625286;
    tfp_guess(13) = 11994.898401;
    tfp_guess(14) = 12356.087390;
    tfp_guess(15) = 12158.123852;
    
end

%% Either load parameters from an existing file, or use the parameters
% listed above
if param_from_file==1
    %     load guess_from_file.txt
    %     guess=guess_from_file;
    fid = fopen('guess_from_file1.txt','r');
    temp = textscan(fid,'%s %s %f');
    fclose(fid);
    guess = temp{3};
else
    if strcmp(cost_specif,'c6') % DOUBLE THRESHOLD
        guess = [alpha; beta; theta; gamma; var_lambda; scale; delta; x_lowbar; x_highbar];
    elseif strcmp(cost_specif,'c5') % LOGISTIC
        guess = [alpha; beta; theta; gamma; var_lambda; scale; delta];
    elseif strcmp(cost_specif,'c7') % EXPONENTIAL
        guess = [beta; theta; gamma; var_lambda; tfp_guess];
    else
        guess = [alpha; beta; theta; gamma; var_lambda; scale; delta];
    end
    
end

%guess(18) = guess(18)*0.95;

disp(['Cost specification: ',cost_specif])

%% CALL THE ROUTINE THAT SOLVES THE MODEL:

if do_estimation == 1
   disp(['Estimation method: ',est_algo]);
    %Setting bounds 
    LB = 0.1*guess; 
    UB = 5.0*guess; 
    if strcmp(cost_specif,'c6')  % DOUBLE THRESHOLD
        %           1      2     3      4        5         6     7       8
        %guess = [alpha; beta; theta; gamma; var_lambda; scale; delta; x_lowbar; x_highbar];
        % keep some parameters fixed at the guess
        LB(4:7) = guess(4:7);
        UB(4:7) = guess(4:7);
        LB(9) = 0.1;
        UB(9) = 0.9;
    elseif strcmp(cost_specif,'c5') % LOGISTIC
        %guess = [alpha; beta; theta; gamma; var_lambda; scale; delta];
        LB(4:7) = guess(4:7);
        UB(4:7) = guess(4:7);
    elseif strcmp(cost_specif,'c1') % QUADRATIC/POLYN
        %guess = [alpha; beta; theta; gamma; var_lambda; scale; delta];
        % estimate only wrt alpha,beta and theta
        LB(1)=0; % alpha
        UB(1)=5;
        LB(3)=1; % theta
        UB(3)=2;
        LB(4:7) = guess(4:7);
        UB(4:7) = guess(4:7);
    elseif strcmp(cost_specif,'c7') % EXPONENTIAL
        %           1      2      3        4        5:19         
        % guess = [beta; theta; gamma; var_lambda; tfp_guess];
        
        %LB(1)=0; UB(1)=0; % beta
        %LB(2)=0; UB(2)=0; % theta
        LB(3)=0.7; UB(3)=0.99; % gamma
        %LB(4)=0; UB(4)=0; % var_lambda
        LB(5:19)=0.8*guess(5:19); UB(5:19)=1.2*guess(5:19); % tfp's
        % plot(year_vec,tfp_guess,year_vec,LB(5:19),'-',year_vec,UB(5:19),'-')
        
    end
    %Bounds are corrected when a guess parameter is negative:
    ind_negative_param = find (guess<0);
    LB(ind_negative_param) = 0.1*guess(ind_negative_param);
    UB(ind_negative_param) = 5.0*guess(ind_negative_param);
    % Check initial guess satisfies bounds
    if any(guess>UB)
        error('guess > UB')
    elseif any(guess<LB)
        error('guess < LB')
    end
    
    
    switch est_algo
        
        case 'nlopt'
            
            % Using NLOPT routines:
            % Setting options
            opt.algorithm = NLOPT_GN_CRS2_LM; %NLOPT_LN_NELDERMEAD;%NLOPT_LN_BOBYQA;%NLOPT_LN_NELDERMEAD;NLOPT_GN_DIRECT_L;%
            opt.lower_bounds = LB;
            opt.upper_bounds = UB;
            opt.min_objective =  @ f_obj;
            % opt.fc = { @Constraint1, @Constraint2 };
            % opt.fc_tol = [1e-5, 1e-5];
            opt.xtol_rel = 1e-5;
            opt.maxeval = 200000;
            opt.maxtime = round(3600*24*6.5);
            %opt.verbose = 1; %Print useful stuff
            % Printing
            fprintf('=========================\n')
            fprintf('ESTIMATION IN PROGRESS...\n')
            fprintf('=========================\n')
            % Call the optimizer
            [optimizer, thevalue, retcode] = nlopt_optimize(opt, guess);
            
        case 'simulan'
            
            %%%%%%%%%%%%%% SIMULATED ANNEALING
            %%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            maxim=0;            % 1 if maximize, 0 otherwise
            rt=0.85;            % temperature decreasing ratio
            eps=1e-2;           % error tolerance for termination
            ns=20;              % number of cycles
            nt=2;
            %max(100,5*npar);   % number of iterations before temperature reduction
            neps=4;             % number of final functions values used to decide upon termination
            maxevl=20000;     % maximum number of function evaluations
            iprint=1;           % control inside printing
            t=1;                % initial temperature
            option=[maxim rt eps ns nt neps maxevl iprint t]';
            c= 2*ones(length(guess),1);
            %             vm = [0.001 0.001 0.001 0.001 0.001 5*ones(1,length(tfp_data))]';
            vm = 0.1*ones(length(guess),1);
            [optimizer,fopt,nacc,nfcnev,nobds,ier,t,vm]= SIMULANS_may(@f_obj,guess,option,LB,UB,c,vm);
            disp('****   results after sa   ****');
            disp('solution');disp(optimizer(:)');
            disp('final step length');disp(vm(:)');
            fprintf('optimal function value              : %g\n',fopt);
            fprintf('number of function evaluations      : %g\n',nfcnev);
            fprintf('number of accepted evaluations      : %g\n',nacc);
            fprintf('number of out-of-bounds evaluations : %g\n',nobds);
            fprintf('final temperature                   : %g\n',t);
            fprintf('error                               : %g\n',ier);
            
        case 'simulan_mathworks'
            
            [x,fval,exitFlag,output] = simulannealbnd(@f_obj,guess,LB,UB);
            fprintf('The number of iterations was : %d\n', output.iterations);
            fprintf('The number of function evaluations was : %d\n', output.funccount);
            fprintf('The best function value found was : %g\n', fval);
            
        case 'fminsearch'
            
            options=optimset('MaxFunEvals',4000,'MaxIter',4000,'FunValCheck','on','TolFun',1e-5,'TolX',1e-5,'disp','on');
            [x,fval,exitflag]=fminsearchcon('f_obj',guess,LB,UB,[],[],[],options);
            
        case 'ga'
            
            %Number of parameter
            n_param = length(guess);
            
            %Objective function
            ObjectiveFunction = @f_obj;
            
            %Constraint function
            %ConstraintFunction = @f_con_ga;
            ConstraintFunction = []; % we don't have any nonlinear constraint!!
            
            %Linear inequality constraints
            Aineq = [];
            Bineq = [];
            
            %Linear equality constraints
            Aeq = [];
            Beq = [];
            
            %Change options
            gaoptions = psoptimset('Cache','on'); %Turn cache on
            gaoptions = gaoptimset(gaoptions,'Display', 'iter'); %Display iterations
            gaoptions = gaoptimset(gaoptions,'TimeLimit',3600*24*3.0); %Time limit of execution
            gaoptions = gaoptimset(gaoptions,'StallTimeLimit',10000); %algorithm stops if no improvement in obj function
            %gaoptio
            ns = gaoptimset(gaoptions,'MutationFcn',@mutationadaptfeasible); %Function for mutation children
            gaoptions = gaoptimset(gaoptions,'FitnessScalingFcn',@fitscalingprop); %Function scales the values of fitness function
            gaoptions=gaoptimset(gaoptions,'PopulationSize',120); %Size of the population
            gaoptions=gaoptimset(gaoptions,'Generations',600); %max number of iterations before the algorithm halts
            if do_parallel==1
                gaoptions = gaoptimset(gaoptions,'UseParallel','always'); %Parallelism
            end
            
            %Printing
            fprintf('=========================\n')
            fprintf('ESTIMATION IN PROGRESS...\n')
            fprintf('=========================\n')
            
            %Call the optimizer
            [x,fval,exitflag,output,population,scores] = ...
                ga(ObjectiveFunction,n_param,Aineq,Bineq,Aeq,Beq,LB,UB,ConstraintFunction,gaoptions);
            
            
    end
   
else
    
    disp('Run the model (NO estimation)')
    %Run at the model:
    if do_experiments==1
        if freeze_te == 1
            disp('T_E kept at 2000 level')
        end
        if freeze_se == 1
            disp('S_E kept at 2000 level')
        end
        if freeze_sw == 1
            disp('S_W kept at 2000 level')
        end
        if freeze_tw == 1
            disp('T_W kept at 2000 level')
        end
        if freeze_productivity == 1
            disp('z(t) kept at 2000 level')
        end
    end
    f_obj(guess);
    
    if do_experiments==1
        experiments
    end
    
    
end




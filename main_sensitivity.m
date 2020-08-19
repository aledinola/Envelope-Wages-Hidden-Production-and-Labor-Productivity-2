%% SENSITIVITY ANALYSIS
% In this script we perform the tax experiments with different values of
% gamma (sensitivity analysis)

%% FOLDER. V11

%% Preliminaries
clear all; clc; close all
format long g

global cost_specif do_estimation year_vec labor_prod_data av_wages informal_size 
global wmin_vec start_year  final_year reform_year
global LB UB est_algo wmin_flag 
global freeze_te freeze_se freeze_sw freeze_tw freeze_productivity experiment_year do_experiments experiment_num


%% Set some flags

do_estimation   = 0;     % if 0, then simply run the model at x = guess  
wmin_flag       = 1;     % set it = 1 if wanna impose the min.wage constraint in the 
                         % NO-tax-evasion problem
                         
%% Set some extra flags for the counterfactual experiments

%-------------- 5 experiments----------------------------------------------

% 1 - keep productivity fixed at 2000 level (==> role of taxes only)
% 2 - keep all taxes fixed
% 3 - change only T_W (other taxes and prod fixed)
% 4 - change only T_E 
% 5 - change only S_E & S_W 
%--------------------------------------------------------------------------

do_experiments = 1; % 1 = do experiments
experiment_num = 1; % we have 5 experiments so far
% If, for ex., freeze_te = 1, then the model is solved keeping the
% corporate business tax "freezed" at its 2000 level.
freeze_te = 1; % corporate business tax
freeze_se = 1; % employers social contribution
freeze_sw = 1; % workers social contribution
freeze_tw = 0; % personal income tax

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

%% We load parameters from an existing file

%     load guess_from_file.txt
%     guess=guess_from_file;
fid = fopen('guess_from_file1.txt','r');
temp = textscan(fid,'%s %s %f');
fclose(fid);
guess = temp{3};

%% We modify "gamma" here
% guess = [beta; theta; gamma; var_lambda; tfp_guess];
% In bechmark, gamma = 0.974102
guess(3) = 0.5;

disp('SENSITIVITY ANALYSIS')

%% CALL THE ROUTINE THAT SOLVES THE MODEL:

disp('Run the model (NO estimation)')
%Run at the model:
if do_experiments==1
    if freeze_te == 1;
        disp('T_E kept at 2000 level')
    end
    if freeze_se == 1;
        disp('S_E kept at 2000 level')
    end
    if freeze_sw == 1;
        disp('S_W kept at 2000 level')
    end
    if freeze_tw == 1;
        disp('T_W kept at 2000 level')
    end
    if freeze_productivity == 1;
        disp('z(t) kept at 2000 level')
    end
end
f_obj(guess);

if do_experiments==1
    experiments
end

    





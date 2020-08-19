%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   "Envelope Wages, Hidden Production and Labor Productivity"            %
%                                                                         %
%       by Alessandro Di Nola, Georgi Kocharkov and A. Vasilev            % 
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all

% This file loads a mat file with counterfactuals, generates figures and tables, 
% and stores them in the subfolder 'results\experiments'

oldpath = fullfile(pwd,'');

% Load results of baseline simulation
load results_baseline

%-------------- 5 experiments----------------------------------------------

% 1 - keep productivity fixed at 2000 level
% 2 - keep all taxes fixed
% 3 - change only T_W (other taxes and prod fixed)
% 4 - change only T_E 
% 5 - change only S_E & S_W 
%--------------------------------------------------------------------------

experiment_num = 1; %choose the experiment

% load counterfactual results
if experiment_num == 1
    load informal_z_fixed 
elseif experiment_num==2
    load informal_alltax_fixed 
elseif experiment_num==3
    load informal_z_fixed 
    load informal_tw_only 
elseif experiment_num==4
    load informal_z_fixed 
    load informal_te_only
elseif experiment_num==5
    load informal_z_fixed 
    load informal_s_only
    
end

% Save plots in subfolder results\experiments

cd 'results\experiments'

if experiment_num == 1
    % EXPERIMENT 1: Keep z(t) fixed at 2000 lovel
    
    plot(year_vec,informal_z_fixed*100,'-o',year_vec,informal_baseline*100,'linewidth',2)
    legend('z(t)=2000','baseline'), grid on
    %title('Informal sector size, over Y')
    xlabel('Year')
    ylabel('Informal Economy as % of GDP')
    print('informal_z_fixed','-depsc')
    print('informal_z_fixed','-dpng')
    savefig('informal_z_fixed')
    
elseif experiment_num==2
    % EXPERIMENT 2: Keep all taxes fixed at 2000 lovel
    
    plot(year_vec,informal_alltax_fixed*100,'-o',year_vec,informal_baseline*100,'linewidth',2)
    legend('taxes=2000','baseline'), grid on
    %title('Informal sector size, over Y')
    xlabel('Year')
    ylabel('Informal Economy as % of GDP')
    print('informal_alltax_fixed','-depsc')
    print('informal_alltax_fixed','-dpng')
    savefig('informal_alltax_fixed')
    
elseif experiment_num==3
    % EXPERIMENT 3: only T_W changes, other taxes are freezed (z changes)
    
    plot(year_vec,informal_tw_only*100,'-o',year_vec,informal_z_fixed*100,'linewidth',2)
    legend('T_{W} only','all taxes','Location','Best'), grid on
    %title('Informal sector size, over Y')
    xlabel('Year')
    ylabel('Informal Economy as % of GDP')
    print('informal_tw_only','-depsc')
    print('informal_tw_only','-dpng')
    savefig('informal_tw_only')
    
elseif experiment_num==4
    % EXPERIMENT 4:  only T_E changes, other taxes are freezed
    
    plot(year_vec,informal_te_only*100,'-o',year_vec,informal_z_fixed*100,'linewidth',2)
    legend('T_{E} only','all taxes'), grid on
    %title('Informal sector size, over Y')
    xlabel('Year')
    ylabel('Informal Economy as % of GDP')
    print('informal_te_only','-depsc')
    print('informal_te_only','-dpng')
    savefig('informal_te_only')
    
    elseif experiment_num==5
    % EXPERIMENT 5:  only S_E and S_W change, other taxes are freezed 
    
    plot(year_vec,informal_s_only*100,'-o',year_vec,informal_z_fixed*100,'linewidth',2)
    legend('S only','all taxes'), grid on
    %title('Informal sector size, over Y')
    xlabel('Year')
    ylabel('Informal Economy as % of GDP')
    print('informal_s_only','-depsc')
    print('informal_s_only','-dpng')
    savefig('informal_s_only')
    
end

%% Table summarizing the findings
cd(oldpath)
load informal_z_fixed 
load informal_alltax_fixed 
cd 'results\experiments'

fid=fopen('results_table.tex','w');
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');

fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Productivity versus Txes - A Decomposition}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lccc}');
fprintf(fid,'\n  \\hline');
fprintf(fid,'\n  & \\textbf{2000} & \\textbf{2014} & \\textbf{Change} \\\\ \\hline');
fprintf(fid,'\n  Data                     &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',36.9, 31,-5.90);
fprintf(fid,'\n  \\textbf{Model}          &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',informal_baseline(1)*100,informal_baseline(end)*100,(informal_baseline(end)-informal_baseline(1))*100);
fprintf(fid,'\n  Change taxes only        &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',informal_z_fixed(1)*100,informal_z_fixed(end)*100,(informal_z_fixed(end)-informal_z_fixed(1))*100);
fprintf(fid,'\n  Change productivity only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',informal_alltax_fixed(1)*100,informal_alltax_fixed(end)*100,(informal_alltax_fixed(end)-informal_alltax_fixed(1))*100);

fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);


%% Go back to the main folder

cd(oldpath)

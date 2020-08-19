%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   "Envelope Wages, Hidden Production and Labor Productivity"            %
%                                                                         %
%       by A. Di Nola, G. Kocharkov and A. Vasilev                        % 
%                                                                         %
%       published on the B.E. Journal of Macroeconomics (2019)            %
%                                                                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all

% This file loads a mat file with results, generates figures and tables, 
% and stores them in the subfolder 'results\model_fit'. The numbering of
% figures and tables refer to the published version of the paper.

load stuff
oldpath = fullfile(pwd,  '' );
cd 'results\model_fit'

%% Generate Figure 10, Section 4.3

% Panel (a): SIZE OF INFORMAL ECONOMY
figure
plot(year_vec,informal_size*100,'-o',year_vec,informal_size_model_reported*100,'linewidth',2)
legend('data','model'), grid on
%title('Informal sector size, over Y hat (reported)')
xlabel('Year')
ylabel('Informal Economy as % of GDP')
print('informal_sector_size','-depsc')
print('informal_sector_size','-dpng')

% Panel (b): OBSERVED WAGES
figure
plot(year_vec,av_wages,'-o',year_vec,w_evasion_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on
ylabel('Real 2005 Euros')
%title('Observed wages')
print('wages','-depsc')
print('wages','-dpng')

% Panel (c): OBSERVED LABOR PRODUCTIVITY
figure
plot(year_vec,labor_prod_data,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on
ylabel('Real 2005 Euros')
%title('Labor Productivity')
print('labor_prod','-depsc')
print('labor_prod','-dpng')

%% Generate Figure 11, Section 4.3

% Panel (a): CROSS-SECTION OF TAX EVASION
figure
plot(1:10,100*h_evasion(:,1)./(tfp_vec(1)*lambda_vec),'-o',...
    1:10,100*h_evasion(:,7)./(tfp_vec(7)*lambda_vec),...
    1:10,100*h_evasion(:,9)./(tfp_vec(9)*lambda_vec),...
    1:10,100*h_evasion(:,15)./(tfp_vec(15)*lambda_vec),'-+','linewidth',2)
axis([-inf inf 0 100])
legend('2000','2006','2008','2014'),grid on
%title('Hidden income: he+hw')
xlabel('Productivity')
print('h_evasion_prod','-depsc')
print('h_evasion_prod','-dpng')

% Panel (b): OBSERVED VS UNOBSERVED LABOR PRODUCTIVITY

figure
plot(year_vec,tfp_vec,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('z(t)','observed','Location','Best'), grid on
%title('Labor Productivity')
print('labor_prod_unobserved','-depsc')
print('labor_prod_unobserved','-dpng')


%% Make table with parameters

%% Generate Table 2, Model Fit, Section 4.3

% Compute average for the period 2000-2014
% Model moments:
av_informal_size_model_reported = mean(informal_size_model_reported);
av_w_evasion_year = mean(w_evasion_year);
av_Y_hat_year = mean(Y_hat_year);

% Data moments:
av_informal_size = mean(informal_size);
av_av_wages = mean(av_wages);
av_labor_prod_data = mean(labor_prod_data);

fid=fopen('model_fit.tex','w');
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');

fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Model Fit - Time-averaged Statistics}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lcc}');
%fprintf(fid,'\n   \\multicolumn{4}{c}{$\\tau_t$} \\\\');
fprintf(fid,'\n  \\hline');
fprintf(fid,'\n  \\textbf{Averages} (2000-2014) & \\textbf{Data} & \\textbf{Model} \\\\ \\hline');
fprintf(fid,'\n  Informal Economy (\\%% of GDP) &  $%5.3f$ & $%5.3f$  \\\\ ',av_informal_size, av_informal_size_model_reported);
fprintf(fid,'\n  Observed Wages              &  $%5.3f$ & $%5.3f$  \\\\ ',av_av_wages, av_w_evasion_year);
fprintf(fid,'\n  Observed Labor Productivity &  $%5.3f$ & $%5.3f$  \\\\ ',av_labor_prod_data, av_Y_hat_year);
fprintf(fid,'\n  Gini Disposable Income      &  $%5.3f$ & $%5.3f$  \\\\ ',0.332, Gini_cw_evasion);
fprintf(fid,'\n  Workers at Min Wage         &  $%5.3f$ & $%5.3f$  \\\\ ',wmin_share_data,wmin_share_model);


fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);

%% Go back to the main folder

cd(oldpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case benchmark, gamma = 0.97
% see main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all
oldpath = fullfile(pwd,  '' );

%% Case benchmark, gamma = 0.97
load gamma097_alltaxes
load gamma097_te_only
load gamma097_tw_only


cd 'results\experiments'
subplot(1,2,1)
plot(year_vec,gamma097_alltaxes,'-o',year_vec,gamma097_te_only,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'FontSize',14)
xlabel('Year','Fontsize',13)

subplot(1,2,2)
plot(year_vec,gamma097_alltaxes,'-o',year_vec,gamma097_tw_only,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)
suptitle('gamma = 0.97 (baseline)')
print('alltax_vs_te_gamma097','-depsc')
print('alltax_vs_te_gamma097','-dpng')


gamma097_alltaxes = gamma097_alltaxes*100;
gamma097_te_only = gamma097_te_only*100;
gamma097_tw_only = gamma097_tw_only*100;

fid=fopen('sensitivity_gamma097_table.tex','w');
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');

fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Decomposing the change in Informality: Benchmark, gamma = 0.97}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lccc}');
%fprintf(fid,'\n   \\multicolumn{4}{c}{$\\tau_t$} \\\\');
fprintf(fid,'\n    & 2000 & 2014 & Change \\\\ \\hline');

fprintf(fid,'\n  All taxes &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma097_alltaxes(1),gamma097_alltaxes(end),gamma097_alltaxes(end)-gamma097_alltaxes(1));
fprintf(fid,'\n  TE only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma097_te_only(1),gamma097_te_only(end),gamma097_te_only(end)-gamma097_te_only(1));
fprintf(fid,'\n  TW only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma097_tw_only(1),gamma097_tw_only(end),gamma097_tw_only(end)-gamma097_tw_only(1));

fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);

%% Case gamma = 0.8
cd(oldpath)
load gamma_alltaxes
load gamma_te_only
load gamma_tw_only

cd 'results\experiments'
subplot(1,2,1)
plot(year_vec,gamma_alltaxes*100,'-o',year_vec,gamma_te_only*100,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'Fontsize',13)
xlabel('Year','Fontsize',13)

subplot(1,2,2)
plot(year_vec,gamma_alltaxes*100,'-o',year_vec,gamma_tw_only*100,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)
suptitle('gamma = 0.8')
print('alltax_vs_te_gamma08','-depsc')
print('alltax_vs_te_gamma08','-dpng')


gamma_alltaxes = gamma_alltaxes*100;
gamma_te_only = gamma_te_only*100;
gamma_tw_only = gamma_tw_only*100;

fid=fopen('sensitivity_gamma08_table.tex','w');
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');

fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Decomposing the change in Informality: Sensitivity Analysis, gamma = 0.8}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lccc}');
%fprintf(fid,'\n   \\multicolumn{4}{c}{$\\tau_t$} \\\\');
fprintf(fid,'\n    & 2000 & 2014 & Change \\\\ \\hline');

fprintf(fid,'\n  All taxes &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma_alltaxes(1),gamma_alltaxes(end),gamma_alltaxes(end)-gamma_alltaxes(1));
fprintf(fid,'\n  TE only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma_te_only(1),gamma_te_only(end),gamma_te_only(end)-gamma_te_only(1));
fprintf(fid,'\n  TW only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma_tw_only(1),gamma_tw_only(end),gamma_tw_only(end)-gamma_tw_only(1));

fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);

%% case gamma = 0.5
cd(oldpath)
load gamma1_alltaxes
load gamma1_te_only
load gamma1_tw_only

cd 'results\experiments'
subplot(1,2,1)
plot(year_vec,gamma1_alltaxes*100,'-o',year_vec,gamma1_te_only*100,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)

subplot(1,2,2)
plot(year_vec,gamma1_alltaxes*100,'-o',year_vec,gamma1_tw_only*100,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',13)
xlabel('Year','Fontsize',13)
suptitle('gamma = 0.5')
print('alltax_vs_te_gamma05','-depsc')
print('alltax_vs_te_gamma05','-dpng')


gamma1_alltaxes = gamma1_alltaxes*100;
gamma1_te_only = gamma1_te_only*100;
gamma1_tw_only = gamma1_tw_only*100;

fid=fopen('sensitivity_gamma05_table.tex','w');
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');

fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Decomposing the change in Informality: Sensitivity Analysis, gamma = 0.5}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lccc}');
%fprintf(fid,'\n   \\multicolumn{4}{c}{$\\tau_t$} \\\\');
fprintf(fid,'\n    & 2000 & 2014 & Change \\\\ \\hline');

fprintf(fid,'\n  All taxes &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma1_alltaxes(1),gamma1_alltaxes(end),gamma1_alltaxes(end)-gamma1_alltaxes(1));
fprintf(fid,'\n  TE only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma1_te_only(1),gamma1_te_only(end),gamma1_te_only(end)-gamma1_te_only(1));
fprintf(fid,'\n  TW only &  $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',gamma1_tw_only(1),gamma1_tw_only(end),gamma1_tw_only(end)-gamma1_tw_only(1));

fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);


%% Yet another figure...

subplot(3,2,1)
plot(year_vec,gamma097_alltaxes,'-o',year_vec,gamma097_te_only,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'FontSize',14)
xlabel('Year','Fontsize',13)

subplot(3,2,2)
plot(year_vec,gamma097_alltaxes,'-o',year_vec,gamma097_tw_only,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)
sgtitle('gamma = 0.97 (baseline)')

subplot(3,2,3)
plot(year_vec,gamma_alltaxes,'-o',year_vec,gamma_te_only,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'Fontsize',13)
xlabel('Year','Fontsize',13)

subplot(3,2,4)
plot(year_vec,gamma_alltaxes,'-o',year_vec,gamma_tw_only,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)
sgtitle('gamma = 0.8')

subplot(3,2,5)
plot(year_vec,gamma1_alltaxes,'-o',year_vec,gamma1_te_only,'linewidth',2)
l=legend('all taxes','T_{E} only'); grid on
set(l,'Fontsize',14)
xlabel('Year','Fontsize',13)

subplot(3,2,6)
plot(year_vec,gamma1_alltaxes,'-o',year_vec,gamma1_tw_only,'linewidth',2)
l=legend('all taxes','T_{W} only'); grid on
set(l,'Fontsize',13)
xlabel('Year','Fontsize',13)
sgtitle('gamma = 0.5')





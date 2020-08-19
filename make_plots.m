%% This script make plots and tables illustrating the MODEL FIT

%save stuff
load stuff
%Save directories
oldpath = fullfile(pwd,  '' );
cd 'results\model_fit'

%% Main figures that we report in the paper (EPS format)
% We don't add the title in Matlab

% INFORMAL SECTOR SIZE
figure
plot(year_vec,informal_size*100,'-o',year_vec,informal_size_model_reported*100,'linewidth',2)
legend('data','model'), grid on
%title('Informal sector size, over Y hat (reported)')
xlabel('Year')
ylabel('Informal Economy as % of GDP')
print('informal_sector_size','-depsc')
print('informal_sector_size','-dpng')

% OBSERVED WAGES
figure
plot(year_vec,av_wages,'-o',year_vec,w_evasion_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on
ylabel('Real 2005 Euros')
%title('Observed wages')
print('wages','-depsc')
print('wages','-dpng')

% OBSERVED LABOR PRODUCTIVITY
figure
plot(year_vec,labor_prod_data,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on
ylabel('Real 2005 Euros')
%title('Labor Productivity')
print('labor_prod','-depsc')
print('labor_prod','-dpng')

% CROSS-SECTION OF TAX EVASION

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

% OBSERVED VS UNOBSERVED LABOR PRODUCTIVITY

figure
plot(year_vec,tfp_vec,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('z(t)','observed','Location','Best'), grid on
%title('Labor Productivity')
print('labor_prod_unobserved','-depsc')
print('labor_prod_unobserved','-dpng')

% MODEL TARGETS TOGETHER

figure
subplot(1,3,1)
plot(year_vec,informal_size*100,'-o',year_vec,informal_size_model_reported*100,'linewidth',2)
legend('data','model'), grid on, axis tight
title('Informal Economy','FontSize',14)
%xlabel('Year')
ylabel('Percentage')
subplot(1,3,2)
plot(year_vec,av_wages,'-o',year_vec,w_evasion_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on, axis tight
title('Wages','FontSize',14)
%xlabel('Year')
ylabel('Real 2005 Euros')
subplot(1,3,3)
plot(year_vec,labor_prod_data,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('data','model','Location','NorthWest'), grid on, axis tight
title('Labor Productivity','FontSize',14)
%xlabel('Year')
ylabel('Real 2005 Euros')
%suptitle('Model Fit')

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'informal_wages_lprod','-dpdf','-r0')

fig = gcf;
fig.PaperPositionMode = 'auto';
print('informal_wages_lprod','-dpng','-r300')
print('informal_wages_lprod','-depsc')


% MODEL GENERATED FACTS

figure
subplot(1,2,1)
plot(1:10,100*h_evasion(:,1)./(tfp_vec(1)*lambda_vec),'-o',...
    1:10,100*h_evasion(:,7)./(tfp_vec(7)*lambda_vec),...
    1:10,100*h_evasion(:,9)./(tfp_vec(9)*lambda_vec),...
    1:10,100*h_evasion(:,15)./(tfp_vec(15)*lambda_vec),'-+','linewidth',2)
axis([-inf inf 0 100])
legend('2000','2006','2008','2014'),grid on
title('Hidden output','FontSize',14)
xlabel('Productivity'),ylabel('Percentage')

subplot(1,2,2)
plot(year_vec,tfp_vec,'-o',year_vec,Y_hat_year,'linewidth',2)
legend('z(t)','observed','Location','Best'), grid on
title('Obs vs Unobs Labor Prod','FontSize',14),ylabel('Real 2005 Euros')

%suptitle('Model-generated Facts')

fig = gcf;
fig.PaperPositionMode = 'auto';
print('modelgen_facts','-dpng','-r0')
print('modelgen_facts','-depsc')


%% Other figures (to inspect some results) png format

figure
plot(year_vec,informal_size,'-o',year_vec,informal_size_model,'linewidth',2)
legend('data','model'), grid on
title('Informal sector size, over Y')
xlabel('Year')
print('informal_sector_size_Y','-dpng')

figure
plot(year_vec,H_year,'-o',year_vec,GDP_year,year_vec,Y_hat_year,'linewidth',2)
legend('Hidden Income','Y','Y hat'), grid on
title('Informal sector size: determinants')
xlabel('Year')
print('informal_sector_size_determinants','-dpng')

figure % declared wages
plot(1:nn_lambda,w_evasion(1:nn_lambda,2),'-o',1:nn_lambda,wmin_vec(2)*ones(nn_lambda,1),'linewidth',2),grid on
legend('w evasion','min wage')

figure % declared profits
plot(1:nn_lambda,e_evasion(1:nn_lambda,tt),'-o','linewidth',2),grid on
legend('e evasion')

% SHARE OF FIRMS AT THE MINIMUM WAGE
figure % share of pairs at the minimum wage
subplot(2,1,1)
plot(year_vec,wmin_year,'-o',year_vec,wmin_share_data*ones(length(year_vec),1),'linewidth',2),grid on
legend('share of pairs s.t. w=wmin (MODEL)','share of pairs s.t. w=wmin (DATA)')
subplot(2,1,2)
plot(year_vec,inactive_year,'-o','linewidth',2),grid on
legend('share of INactive pairs')
print('share_minwage','-dpng')

figure
plot(1:10,he_evasion(:,1)./(tfp_vec(1)*lambda_vec),'-o',...
    1:10,he_evasion(:,7)./(tfp_vec(7)*lambda_vec),...
    1:10,he_evasion(:,9)./(tfp_vec(9)*lambda_vec),...
    1:10,he_evasion(:,15)./(tfp_vec(15)*lambda_vec),'-+','linewidth',2)
legend('2000','2006','2008','2014'),grid on
title('Hidden income he')
xlabel('Productivity')
print('he_polfun','-dpng')

figure
plot(1:10,hw_evasion(:,1)./(tfp_vec(1)*lambda_vec),'-o',...
    1:10,hw_evasion(:,7)./(tfp_vec(7)*lambda_vec),...
    1:10,hw_evasion(:,9)./(tfp_vec(9)*lambda_vec),...
    1:10,hw_evasion(:,15)./(tfp_vec(15)*lambda_vec),'-+','linewidth',2)
legend('2000','2006','2008','2014'),grid on
title('Hidden income hw')
xlabel('Productivity')
print('hw_polfun','-dpng')

figure
plot(year_vec,y_evasion_year,year_vec,Y_hat_year)
legend('decisions','resource constraint')
title('Labor Productivity: CHECK')
% note: if E(lambda)=1, the two measures must coincide

%% Make table with parameters

%% Make table with aggregate results

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
fprintf(fid,'\n \\caption{Aggregate Statistics}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\label{results1}');
fprintf(fid,'\n \\begin{tabular}{lcc}');
%fprintf(fid,'\n   \\multicolumn{4}{c}{$\\tau_t$} \\\\');
fprintf(fid,'\n   Averages (2000-2014) & Data & Model \\\\ \\hline');
fprintf(fid,'\n  Informal sector &  $%5.3f$ & $%5.3f$  \\\\ ',av_informal_size, av_informal_size_model_reported);
fprintf(fid,'\n  Observed wages &  $%5.3f$ & $%5.3f$  \\\\ ',av_av_wages, av_w_evasion_year);
fprintf(fid,'\n  Output per worker &  $%5.3f$ & $%5.3f$  \\\\ ',av_labor_prod_data, av_Y_hat_year);
fprintf(fid,'\n  Gini &  $%5.3f$ & $%5.3f$  \\\\ ',0.332, Gini_cw_evasion);
fprintf(fid,'\n  share of establ at minwage &  $%5.3f$ & $%5.3f$  \\\\ ',wmin_share_data,wmin_share_model);


fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

fprintf(fid,'\n \\end{document}');
fclose(fid);

%% Go back to the main folder

cd(oldpath)

%% save results for later use
% legend: varname_baseline

informal_baseline = informal_size_model_reported;
save results_baseline informal_baseline year_vec



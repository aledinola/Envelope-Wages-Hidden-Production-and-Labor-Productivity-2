function [ tax ] = T_W( wage )
%This function computes the labor income tax payments of the worker

global time display_results
global T_W_year cubic interpol
global freeze_tw experiment_year do_experiments

if do_experiments==1 && freeze_tw==1
    year = experiment_year;
else
    year = time;
end


if interpol==1
   
    indt = year-2000+1; % 2000 is 1, 2001 is 2, ..., 2014 is 15
    T_W_year_t = T_W_year{indt}; % 2-column vector
    if display_results==1
        % check that wage falls in the range of T_W_year_t(:,1)
        if wage < T_W_year_t(1,1)  
            %error('wage out of range!')
            disp('LABOR INCOME TAX: wage < min wage!')
        elseif wage > T_W_year_t(end,1)
            disp('LABOR INCOME TAX: wage > max wage')
        end
    end
    
    if cubic==1
        tax = max(10^-5,pchip(T_W_year_t(:,1),T_W_year_t(:,2),wage));
    else
        % using Matlab predefined routine for linear interpolation
        tax = interp1(T_W_year_t(:,1),T_W_year_t(:,2),wage,'linear','extrap');
                      % X        % Y         % xi
     end                 
else
    
    %*********************OLD CODE *****************************************
    switch year
    
    case 2000
        tax = ( 0.117 + 0.088 * log(wage/wagebar) ) * wage;
    case 2001
        tax = ( 0.107 + 0.096 * log(wage/wagebar) ) * wage;
    case 2002
        tax = ( 0.098 + 0.087 * log(wage/wagebar) ) * wage;
    case 2003
        tax = ( 0.096 + 0.083 * log(wage/wagebar) ) * wage;
    case 2004
        tax = ( 0.096 + 0.082 * log(wage/wagebar) ) * wage;
    case 2005
        tax = ( 0.087 + 0.065 * log(wage/wagebar) ) * wage;
    case 2006
        tax = ( 0.074 + 0.073 * log(wage/wagebar) ) * wage;
    case 2007
        tax = ( 0.081 + 0.076 * log(wage/wagebar) ) * wage;
    case 2008
        tax = 0.100 * wage; 
    case 2009
        tax = 0.100 * wage; 
    case 2010
        tax = 0.100 * wage; 
    case 2011
        tax = 0.100 * wage; 
    case 2012
        tax = 0.100 * wage; 
    case 2013
        tax = 0.100 * wage; 
    case 2014
        tax = 0.100 * wage; 
        
    end

end

end

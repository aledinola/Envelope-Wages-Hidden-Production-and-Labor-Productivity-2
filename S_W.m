function [ tax ] = S_W( wage )
%This function computes the social security payments of the worker

global time display_results

global S_W_year cubic interpol
global freeze_sw experiment_year do_experiments


if do_experiments==1 && freeze_sw==1
    year = experiment_year;
else
    year = time;
end

if interpol==1
    
        indt = year-2000+1; % 2000 is 1, 2001 is 2, ..., 2014 is 15
        S_W_year_t = S_W_year{indt}; % 2-column vector
        if display_results==1
            % check that wage falls in the range of S_W_year_t(:,1)
            if wage < S_W_year_t(1,1)  
                %error('wage out of range!')
                disp('wage < min wage!')
            elseif wage > S_W_year_t(end,1)
                disp('wage > max wage')
            end
        end
        
        if cubic==1
            tax = max(10^-5,pchip(S_W_year_t(:,1),S_W_year_t(:,2),wage));
        else
            % using Matlab predefined routine for linear interpolation
            tax = interp1(S_W_year_t(:,1),S_W_year_t(:,2),wage,'linear','extrap');
                          % X        % Y         % xi
        end
        
else

        switch year

            case 2000
                tax = 0.1205 * wage;
            case 2001
                tax = 0.1130 * wage;
            case 2002
                tax = 0.1080 * wage;
            case 2003
                tax = 0.1155 * wage;
            case 2004
                tax = 0.1080 * wage;
            case 2005
                tax = 0.1245 * wage;
            case 2006
                tax = 0.1243 * wage;
            case 2007
                tax = 0.1243 * wage;
            case 2008
                tax = 0.1300 * wage;
            case 2009
                tax = 0.1300 * wage;
            case 2010
                tax = 0.1210 * wage;
            case 2011
                tax = 0.1290 * wage;
            case 2012
                tax = 0.1290 * wage;
            case 2013
                tax = 0.1290 * wage;
            case 2014
                tax = 0.1290 * wage;

        end

end

end


function [ tax ] = S_E( wage )
%This function computes the social security payments of the employer
% Two methods:
% (1) - interpolation using the true tax data
% (2) - estimated log specification (see Guner et al.)
global time display_results
global S_E_year cubic interpol
global freeze_se experiment_year do_experiments


if do_experiments==1 && freeze_se==1
    year = experiment_year;
else
    year = time;
end


if interpol==1
 
        
   
    indt = year-2000+1; % 2000 is 1, 2001 is 2, ..., 2014 is 15
        S_E_year_t = S_E_year{indt}; % 2-column vector

        % plot(S_E_year_t(:,1),S_E_year_t(:,2),'o')
        if display_results==1
            % check that wage falls in the range of S_E_year_t(:,1)
            if wage < S_E_year_t(1,1) 
                %error('wage out of range!')
                disp('wage < min wage!')
            elseif wage > S_E_year_t(end,1)
                 disp('wage > max wage!')
            end
        end

        if cubic==1
            %Piecewise Cubic Hermite Interpolating Polynomial 
            tax = max(10^-5,pchip(S_E_year_t(:,1),S_E_year_t(:,2),wage));
        else
            % linear
            tax = interp1(S_E_year_t(:,1),S_E_year_t(:,2),wage,'linear','extrap');
        end

else
        switch year

            case 2000
                tax = 0.3615 * wage;
            case 2001
                tax = 0.3390 * wage;
            case 2002
                tax = 0.3465 * wage;
            case 2003
                tax = 0.3165 * wage;
            case 2004
                tax = 0.3110 * wage;
            case 2005
                tax = 0.2410 * wage;
            case 2006
                tax = 0.2332 * wage;
            case 2007
                tax = 0.2247 * wage;
            case 2008
                tax = 0.1740 * wage;
            case 2009
                tax = 0.1690 * wage;
            case 2010
                tax = 0.1590 * wage;
            case 2011
                tax = 0.1580 * wage;
            case 2012
                tax = 0.1580 * wage;
            case 2013
                tax = 0.1580 * wage;
            case 2014
                tax = 0.1580 * wage;

        end
end
end


function [ tax ] = T_E( profit)
%This function computes the corpotate tax payments of the employer


global time
global freeze_te experiment_year do_experiments

if do_experiments==1 && freeze_te==1
    year = experiment_year;
else
    year = time;
end
 
switch year
    
    case 2000
        tax = 0.325 * profit;
    case 2001
        tax = 0.280 * profit;
    case 2002
        tax = 0.235 * profit;
    case 2003
        tax = 0.235 * profit;
    case 2004
        tax = 0.195 * profit;
    case 2005
        tax = 0.150 * profit;
    case 2006
        tax = 0.150 * profit;
    case 2007
        tax = 0.100 * profit;
    case 2008
        tax = 0.100 * profit;
    case 2009
        tax = 0.100 * profit;
    case 2010
        tax = 0.100 * profit;
    case 2011
        tax = 0.100 * profit;
    case 2012
        tax = 0.100 * profit;
    case 2013
        tax = 0.100 * profit;
    case 2014
        tax = 0.100 * profit;
        
end

end





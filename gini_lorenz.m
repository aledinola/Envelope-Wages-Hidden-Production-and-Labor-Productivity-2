function [F,L,gini] = gini_lorenz(RelInc,Dist)

% This function computes Lorenz curve and Gini coefficient for a population with a discrete
% density function
% Lorenz curve with full equality line: plot(F,F,F,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -you need a vector of incomes, RelInc
% -you need a vector of probability mass for each income, Dist

%Sort into ascending order
%Rank the relative incomes and derive their indices
[r_ascend, r_ascend_index]= sort(RelInc,'ascend');

%Form the deciles in terms of cutoffs in the vectors f and r
f_ascend = Dist(r_ascend_index); %ordered f's by ascending r's
F_cum = cumsum(f_ascend); %cumulative ordered f's by ascending r's
C_cum = cumsum(f_ascend.*r_ascend); %cumulative share of income
%Calculate deciles
% perc10 = wprctile(r_ascend,10,f_ascend);
% perc50 = wprctile(r_ascend,50,f_ascend);
% perc90 = wprctile(r_ascend,90,f_ascend);
% 
% %Ratios
% ratio_90_10 = perc90 / perc10;
% ratio_90_50 = perc90 / perc50;
% ratio_50_10 = perc50 / perc10;

%Calculate Gini   WRONG
gini = 0;
for j=1: size(F_cum,1)-1
   
    gini =gini + F_cum(j)*C_cum(j+1) - F_cum(j+1)*C_cum(j);
   
end

F = F_cum;
L = C_cum/C_cum(end);

%%%%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = length(y);
% n1 = length(prob);
% 
% if n~=n1
%     error('Check GINI LORENZ')
% end
% 
% F1 = cumsum(prob); % column vector n*1
% F = [0;F1];
% S = zeros(n,1);
% S(1) = y(1)*prob(1);
% for i=2:n
%     S(i) = S(i-1)+y(i)*prob(i);
% end
% L1 = zeros(n,1);
% for i=1:n
%     L1(i) = S(i)/S(n);
% end
% L = [0;L1];
% 
% % plot(F,F,F,L)

end


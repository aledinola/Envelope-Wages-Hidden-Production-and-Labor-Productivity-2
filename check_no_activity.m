%This m-file computes the region of productivity where pairs are inactive.

% When drop_points = [], i.e. no points in lambda_vec should be
% dropped

% When drop_points is not empty, it contains the indices of lambda_vec to
% be dropped

%Step 1: check if the lowest productivity point can afford to pay the
%minimum wage

if tfp * lambda_vec(1) >= wmin 
       
    drop_points = [];
    
else
    
    %Step 2; find j_bar, i.e. the highest productivity point where they
    %cannot pay the minimum wage
    
    j_bar = find(wmin<tfp*lambda_vec,1) - 1;
    
    %Step 3: find i_bar, i.e. this is the first (the lowest) point in 
    %lambda_vec at which the cost of hiding everything cannot be covered 
    
    i_bar = find( f_kappa(tfp*lambda_vec,0)<=tfp*lambda_vec, 1,'last') + 1;
    
    %Step 4: find the inactivity region
    drop_points  = i_bar:j_bar;
    
end
    
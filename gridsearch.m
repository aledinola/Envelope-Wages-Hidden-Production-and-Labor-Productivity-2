% Computes "optimal" (he,hw) on a coarse grid and 
% gives as output a refined guess
L = n_grid_init;
 
h_vec = linspace(0,(tfp*lambda),L)';
% if t>=2
%     h_past = h_evasion(n,t-1);
%     hopt = h_evasion_prod(n,t-1)*tfp*lambda;
%     index=max(find(hopt>=h_vec));
%     h_vec = [h_vec(1:index); hopt; h_vec(index+1:L-1)];
%     %index1 = max(find(h_past>=h_vec));
%     %h_vec = [h_vec(1:index1); h_past; h_vec(index1+1:L-1)];  
% end
obj_evasion_vec = zeros(L,1);

for i=1:L
    hh = h_vec(i);
    obj_evasion_vec(i) = obj_evasion(hh,param_cost);
    % if constraints are violated, than obj = INF
    % taken care inside "obj_evasion" function
end

if length(find(obj_evasion_vec==inf))==L
    warning('There are no feasible points to be chosen as initial guess!')
    keyboard
end


[maxobj,h_ind_guess] = min(obj_evasion_vec);
guess = h_vec(h_ind_guess);
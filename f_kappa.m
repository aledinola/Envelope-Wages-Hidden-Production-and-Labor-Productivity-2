function [ cost ] = f_kappa(h,param_cost)
%This function computes the cost of hiding income

% global alpha delta x_bar x_lowbar x_highbar
global time cost_specif
% global lambda_vec h_evasion t n cost_specif

%alpha = param_cost.alpha; 
beta = param_cost.beta;
theta = param_cost.theta;

switch cost_specif
    
    case 'c1' % alpha+beta*h^theta + delta(h/y-h_past/y_past)^2
        if time>2000
            y = tfp_vec(t)*lambda_vec(n);
            y_past = tfp_vec(t-1)*lambda_vec(n);
            h_past = h_evasion(n,t-1);
            cost  = alpha + beta * (h).^theta + delta*((h/y)-(h_past/y_past)).^2;
            %cost  = delta*((h/y)-(h_past/y_past)).^2; % only adj costs, as
            %a check
        else
            cost  = alpha + beta * (h).^theta;
        end
    case 'c2'
        if time>2000
            y = tfp_vec(t)*lambda_vec(n);
            y_past = tfp_vec(t-1)*lambda_vec(n);
            h_past = h_evasion(n,t-1);
            cost  = alpha + beta * (h).^theta + delta*((h/y)-(h_past/y_past)).^2;
            %cost  = delta*((h/y)-(h_past/y_past)).^2; % only adj costs, as
            %a check
        else
            cost  = alpha + beta * (h).^theta;
        end
        
    case 'c3' % cost depends on h/y ==> high prod evade more
        y = tfp_vec(t)*lambda_vec(n);
        
        if time>2000
            y_past = tfp_vec(t-1)*lambda_vec(n);
            h_past = h_evasion(n,t-1);
            cost  = alpha + beta *(h/y).^theta + delta*((h/y)-(h_past/y_past)).^2;
        else
            cost  = alpha + beta * y*(h).^theta;
        end
        
    case 'c4' % beta * (alpha*h)^theta/(tfp_vec(t)^theta)
        y = tfp_vec(t)*lambda_vec(n);
        
        if time>2000
            y_past = tfp_vec(t-1)*lambda_vec(n);
            h_past = h_evasion(n,t-1);
            cost  = beta * (alpha*h)^theta/(tfp_vec(t)^theta)+ delta*((h/y)-(h_past/y_past)).^2;
        else
            cost  = beta * (alpha*h)^theta/(tfp_vec(t)^theta);
        end
        
    case 'c5' % LOGISTIC alpha/1+exp[-theta(h-beta)]
        
        if time>2000
            %y_past = tfp_vec(t-1)*lambda_vec(n);
            h_past = h_evasion(n,t-1);
            cost  = alpha ./ (1+exp(-theta*(h-beta))) + delta*((h)-(h_past)).^2;
        else
            cost  = alpha ./ (1+exp(-theta*(h-beta)));
        end
        
    case 'c6' % double threshold
        
        y = tfp_vec(t)*lambda_vec(n);
        
        if h/y > x_highbar
            cost = alpha+beta*y^theta*(x_highbar-x_lowbar)^theta;
        elseif h/y < x_lowbar
            cost = 0;
        else
            cost  = alpha + beta *(h-x_lowbar*y).^theta;
        end
    case 'c7' % exponential
        %cost = alpha+beta*exp(theta*(h-x_bar));
       
        cost = beta*exp(theta*h);
        
end




end


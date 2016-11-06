function [H,grad,iter] = nlssubprob_accel(V,W,Hinit,tol,maxiter,W_constraints)
% H, grad: output solution and gradient
% iter: #iterations used
% Input:
%   V, W: constant matrices
%   Hinit: initial solution
%   tol: stopping tolerance
%   maxiter: limit of iterations
%   W_constraints: type of constraint on W {'positive','on_simplex','inside_simplex'}
%

H = Hinit; 
WtV = W'*V;
WtW = W'*W; 

    L=norm(WtW);	% Lipschitz constant
%     H=Z;    % Initialization
    H_accel=H;
    grad = WtW*H - WtV;
    alpha1=1;

    for iter=1:maxiter,
        H0=H;
        H = project_proportions( H_accel - grad/L, W_constraints);
%         H=max(H_accel-grad/L,0);    % Calculate sequence 'Y'
        alpha2=0.5*(1+sqrt(1+4*alpha1^2));
        H_accel=H + ((alpha1-1)/alpha2)*(H-H0);
        alpha1=alpha2;
        grad = WtW*H - WtV;
        
      stop_crit = get_stopping_cond(H_accel,grad, W_constraints);
      if stop_crit < tol,
        break
      end
    end    
    

if iter==maxiter,
   fprintf('Max iter in nlssubprob\n');
end
end


function H = project_proportions( H, projection_type)
% project the proportion onto a set
%   'on_simplex' - projects the W such that each row sums to one
%   'inside_simplex' - projects W such that each row sums is in the interval [0,1]
%   'positive' - make sure each value in W is non-negative
%   'on_simplex_with_noise' - projects the W such that each row sums to one
%        and adds an eps vector such that there would not be a vector which
%        is all zeros

    switch projection_type
        case 'on_simplex'
            H = stochasticMatrixProjection(H,'on');
        case 'on_simplex_with_noise'
            H = stochasticMatrixProjection(H,'on');
            H = H + rand(size(H))* eps; % moves it away from zero
        case 'inside_simplex'
            H = stochasticMatrixProjection(H,'inside');
        case 'positive'
            H=(H>0).*H;
        otherwise
            error('unknown option for projection type - %s', projection_type);
    end

end

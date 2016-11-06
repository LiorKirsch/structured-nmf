function [H,grad,iter] = nlssubprob(V,W,Hinit,tol,maxiter,W_constraints,accelerated)
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations
% accelerated: do nesterov accelerated gradient

H = Hinit; 
WtV = W'*V;
WtW = W'*W; 

alpha = 1; beta = 0.1;

if accelerated
    
    
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
    
else
    for iter=1:maxiter,  
      grad = WtW*H - WtV;

      stop_crit = get_stopping_cond(H,grad, W_constraints);
      if stop_crit < tol,
        break
      end

      % projected gradient with step-size search 
      [H, alpha] = proj_grad_with_s_search(H,WtW,grad,alpha, beta,W_constraints);
    end
end
if iter==maxiter,
   fprintf('Max iter in nlssubprob\n');
end
end

function [H, alpha] = proj_grad_with_s_search(H,WtW,grad,alpha, beta,W_constraints)

 Hn = project_proportions( H - alpha*grad, W_constraints);
 
  d = Hn-H;
  gradd=sum(sum(grad.*d));  % grad(f(x_k)) * (x_k+1 - x_k)
  dQd = sum(sum((WtW*d).*d));
  if gradd + 0.5*dQd > 0.01*gradd, 
    % decrease alpha
    while 1,
      alpha = alpha*beta;
      
      Hn = project_proportions( H - alpha*grad, W_constraints);

      d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd <= 0.01*gradd | alpha < 1e-20,      
        H = Hn; break;
      end
    end 
  else 
    % increase alpha
    while 1,
      Hp = Hn;
      alpha = alpha/beta;
      
      Hn = project_proportions( H - alpha*grad, W_constraints);
      
      d = Hn-H;
      gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
      if gradd + 0.5*dQd > 0.01*gradd | Hn == Hp | alpha > 1e10,      
        H = Hp; alpha = alpha*beta; break;
      end
    end 
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


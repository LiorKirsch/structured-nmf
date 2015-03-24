function do_early_stop = check_for_early_stopping(X,W,H,init_time,initgrad,parms)
    do_early_stop = false;
    tol = take_from_struct(parms, 'tolerance', 10^-6);
    timelimit = take_from_struct(parms, 'timelimit', 10000);

    gradW = W*(H*H') - X*H';
    gradH = (W'*W)*H - W'*X;
    
    % KKT conditions:
    % gradW >= 0
    % W >= 0
    % W *. gradW = 0;   ===== same as gradW(W>0) = 0
    projnorm = norm( [norm(gradW(gradW<0 | W>0.00001)); norm(gradH(gradH<0 | H>0.00001))] );
    if (projnorm<tol*initgrad || (cputime-init_time) > timelimit)
      fprintf('Early stoping Final proj-grad norm %f', projnorm);
      do_early_stop = true;
    end

end

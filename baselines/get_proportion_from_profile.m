function W = get_proportion_from_profile(X,H, parms)

    als_solver = 'pinv_project';

    % Minimzation step.
    switch als_solver
        case 'pinv_project'
            W = ((pinv(H*H')*H)*X')';
%===== I removed this because this is done in a next when I project to or on the simplex.
%             W=(W>0).*W;  
%=====
        case 'active_set'
            [W,gradW,subIterW] = nnlsm_activeset(H',X',1,0,init_W'); 
            W=W'; 
        case 'blockpivot'
            [W,gradW,subIterW] = nnlsm_blockpivot(H',X',0,init_W');
            W=W'; 
    end
    
    % Projection step.
    switch parms.W_constrains
        case 'on_simplex'
            W = stochasticMatrixProjection(W');
            W = W';
        case 'inside_simplex'
            W=(W>0).*W;
            row_sum = sum(W,2);
            valid_rows =  (0 <= row_sum) & (row_sum <= 1);
            W_tmp = stochasticMatrixProjection(W(~valid_rows,:)');
            W(~valid_rows,:) = W_tmp';
        case 'positive'
             if strcmp(als_solver, 'pinv_project')
                W=(W>0).*W;
            end
        otherwise
            error('unknown option for w_constrains - %s', W_constrains);
    end
end
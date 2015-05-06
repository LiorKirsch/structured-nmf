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
    W = project_proportions( W, W_constraints );
    
end
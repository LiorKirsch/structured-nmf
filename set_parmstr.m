function parmstr = set_parmstr(parms)
%
%
%

    parmstr = '';
    
    if isfield(parms, 'num_types')
        parmstr = sprintf('%s_K%d', parmstr, parms.num_types);
    end    
    
    if isfield(parms, 'num_samples')
        parmstr = sprintf('%s_N%d', parmstr, parms.num_samples);
    end
    
    if isfield(parms, 'maxiter')
        parmstr = sprintf('%s_IT%d', parmstr, parms.maxiter);
    end
        
    if isfield(parms, 'W_constrains')
        parmstr = sprintf('%s_%s', parmstr, parms.W_constrains);
    end    
    
    if isfield(parms, 'rand_seed')
        parmstr = sprintf('%s_S%d', parmstr, parms.rand_seed);
    end        
    
    if isfield(parms, 'num_restarts')
        parmstr = sprintf('%s_NN%d', parmstr, parms.num_restarts);
    end    

    if isfield(parms, 'H_lambda')
        parmstr = sprintf('%s_HL%g', parmstr, parms.H_lambda);
    end    

    if isfield(parms, 'W_lambda')
        parmstr = sprintf('%s_WL%g', parmstr, parms.W_lambda);
    end    
    
    if isfield(parms, 'subsample_iter')
        parmstr = sprintf('%s_SSI%d', parmstr, parms.subsample_iter);
    end    

    if isfield(parms, 'subsample_repeats')
        if parms.subsample_repeats ~= 10
            parmstr = sprintf('%s_SSR%d', parmstr, parms.subsample_repeats);
        end
    end    
    
    if isfield(parms, 'nmf_method')
        parmstr = sprintf('%s_%s', parmstr, parms.nmf_method);
    end 

    if isfield(parms, 'log_transform')
        if parms.log_transform
            parmstr = sprintf('%s_log', parmstr);
        end
    end 
    

end
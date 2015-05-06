function [parmstr, dirparmstr] = set_parmstr(parms)
%
%
%

    parmstr = '';
    dirparmstr = '';
    
    if isfield(parms, 'num_types')
        parmstr = sprintf('%s_K%d', parmstr, parms.num_types);
    end    
    
    if isfield(parms, 'num_samples')
        parmstr = sprintf('%s_N%d', parmstr, parms.num_samples);
    end
    
    if isfield(parms, 'maxiter')
        parmstr = sprintf('%s_IT%d', parmstr, parms.maxiter);
    end
        
    if isfield(parms, 'W_constraints')
        parmstr = sprintf('%s_%s', parmstr, parms.W_constraints);
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

    if isfield(parms, 'multi_region')
        if ~strcmp(parms.multi_region,'');
            parmstr = sprintf('%s_%s', parmstr, parms.multi_region);
        end
    end 
    
    if isfield(parms, 'prior_dataset')
        parmstr = sprintf('%s_Pri_%s', parmstr,parms.prior_dataset);
        parmstr = sprintf('%s_%s', parmstr,strjoin(parms.prior_types,'_') );
    end 
    
    if isfield(parms, 'structre_type')
        parmstr = sprintf('%s_%s', parmstr, parms.structre_type);
    end 


    if isfield(parms, 'mix_files')
        mix_files = sort(parms.mix_files);
        dirparmstr = fullfile(dirparmstr,strjoin(mix_files,'/'));
    end
    if isfield(parms, 'dataset_file')
        dirparmstr = parms.dataset_file;
    end
    if isfield(parms, 'dataset_file') && isfield(parms, 'mix_files')
        error('cannot specifiy both "dataset_file" and "mix_files"');
    end
end
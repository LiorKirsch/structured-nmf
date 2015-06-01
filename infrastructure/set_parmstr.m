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
    
    if isfield(parms, 'init_type')
        if strcmp(parms.init_type, 'random')
            if isfield(parms, 'rand_seed')
                parmstr = sprintf('%s_S%d', parmstr, parms.rand_seed);
            end
            if isfield(parms, 'num_restarts')
                parmstr = sprintf('%s_NN%d', parmstr, parms.num_restarts);
            end    
        else
            parmstr = sprintf('%s_init%s', parmstr, parms.init_type);
        end
    else
        if isfield(parms, 'rand_seed')
            parmstr = sprintf('%s_S%d', parmstr, parms.rand_seed);
        end
        if isfield(parms, 'num_restarts')
            parmstr = sprintf('%s_NN%d', parmstr, parms.num_restarts);
        end    
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
        
        switch parms.structre_type
          case 'tree'

          case {'relations','relations_dist','relations_parentdist',...
                'relations_parent_level'}
            %add a hash for the matrix
            sum_hash = sum(sum(parms.structure_matrix));
            parmstr = sprintf('%s_Hash%g', parmstr, sum_hash);
            
            parms.
            if isfield(parms, 'structure_filter')
                parmstr = sprintf('%s_fltr%g', parmstr, parms.structure_filter);
            end
            
            if isfield(parms, 'do_sep_init')
                if parms.do_sep_init
                    parmstr = sprintf('%s_sepInit', parmstr);
                end
            end
            otherwise
                error('unkown strcuture type %s',parms.structre_type);
        end
    end 

    if isfield(parms,'num_markers')
        parmstr = sprintf('%s_nummrk%d', parmstr, parms.num_markers);
    end
    
    if isfield(parms, 'H_markers')
        if iscell(parms.H_markers)
            sum_hash = 0;
            for i_cell = 1:length(parms.H_markers)
                tmp_matrix = reshape(1:numel(parms.H_markers{i_cell}), size(parms.H_markers{i_cell}));
                tmp_matrix = tmp_matrix .* parms.H_markers{i_cell};
                sum_hash = sum_hash +    sum(tmp_matrix(:)) ;
            end
        else
            tmp_matrix = reshape(1:numel(parms.H_markers), size(parms.H_markers));
            tmp_matrix = tmp_matrix .* parms.H_markers;
            sum_hash = sum(tmp_matrix(:)) ;
        end
        parmstr = sprintf('%s_MrkHash%g', parmstr, sum_hash);
    end 

    if isfield(parms, 'gene_subset')
        if ~strcmp( parms.gene_subset,'all')
            parmstr = sprintf('%s_%s', parmstr, parms.gene_subset);
            if isfield(parms, 'gene_hash')
                parmstr = sprintf('%s_GenesHash%g', parmstr, parms.gene_hash);
            end
            if isfield(parms, 'gene_okaty_filter')
                parmstr = sprintf('%s_fltr%s', parmstr, parms.gene_okaty_filter);
            end
        end
        
    end
    
    if isfield(parms, 'H_lambda_priors')
        lambda_prior_string = sprintf('%g_', parms.H_lambda_priors);
        parmstr = sprintf('%s_lambdaPrr%s', parmstr, lambda_prior_string);
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
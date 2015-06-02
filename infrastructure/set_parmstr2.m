function [parmstr, dirparmstr] = set_parmstr2(parms)
%    
    parmstr = ''; dirparmstr = '';
        
    parmstr = add_parameter('num_types', 'K%d', parms, parmstr);
    parmstr = add_parameter('num_samples', 'N%d', parms, parmstr);
    parmstr = add_parameter('maxiter', 'IT%d', parms, parmstr);
    parmstr = add_parameter('W_constraints', '%s', parms, parmstr);
    
    if isfield(parms, 'init_type')
        if strcmp(parms.init_type, 'random')
            parmstr = add_parameter('rand_seed', 'S%d', parms, parmstr);
            parmstr = add_parameter('num_restarts', 'NN%d', parms, parmstr);
        else
            parmstr = add_parameter('init_type', 'init%s', parms, parmstr);
        end
    else
        parmstr = add_parameter('rand_seed', 'S%d', parms, parmstr);
        parmstr = add_parameter('num_restarts', 'NN%d', parms, parmstr);        
    end        
    
    parmstr = add_parameter('H_lambda', 'HL%g', parms, parmstr);    
    parmstr = add_parameter('W_lambda', 'WL%g', parms, parmstr);
    parmstr = add_parameter('subsample_iter', 'SSI%d', parms, parmstr);
    parmstr = add_parameter('subsample_repeats', 'SSR%d', parms, parmstr, 10);        
    parmstr = add_parameter('nmf_method', '%s', parms, parmstr);    

    if isfield(parms, 'log_transform')
        if parms.log_transform
            parmstr = sprintf('%s_log', parmstr);
        end
    end
    parmstr = add_parameter('multi_region', '%s', parms, parmstr, '');
    
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
            % Add a hash for the matrix
            sum_hash = sum(sum(parms.structure_matrix));
            parmstr = sprintf('%s_Hash%g', parmstr, sum_hash); ...
                      
  
            if isfield(parms, 'structure_filter')
                if ~isnan(parms.structure_filter)
                    parmstr = sprintf('%s_fltr%g', parmstr, parms.structure_filter);
                end
            end                                
            parmstr = add_parameter('do_sep_init', 'sepInit', parms, parmstr, false);
          otherwise
            error('unkown strcuture type %s',parms.structre_type);
        end
    end 
    parmstr = add_parameter('num_markers', 'nummrk%d', parms, parmstr);
    
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
        if ~strcmp(parms.gene_subset,'all')
            parmstr = sprintf('%s_%s', parmstr, parms.gene_subset);
            parmstr = add_parameter('gene_hash', 'GenesHash%g', parms, parmstr);
            parmstr = add_parameter('gene_okaty_filter', 'fltr%s', parms, parmstr);
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
    
    
    parmstr = add_parameter('random_init_spread', 'RIS%g', parms, parmstr, 1)
    
end

% ========================================================
function parmstr = add_parameter(fieldname, format, ...
                                 parms, parmstr, value_to_skip)
    if ~isfield(parms, fieldname)
        return
    end    
    if exist('value_to_skip', 'var') && parms.(fieldname) == value_to_skip
        return
    end
    if findstr('%', format)
        parmstr = sprintf(['%s_' format], parmstr, parms.(fieldname));
    else        
        parmstr = [parmstr '_' format];
    end
end

function parms = conf(parms)
%
% 'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise';
% alg_list = {'alsPinv', 'alsActiveSet', 'mm', 'alsBlockpivot','cjlin', 'prob'};     
%
    [~, parms] = take_from_struct(parms, 'num_types', 3);  
    [~, parms] = take_from_struct(parms, 'maxiter', 500);
    [~, parms] = take_from_struct(parms, 'loglevel', 0);  
    [~, parms] = take_from_struct(parms, 'W_constraints', 'on_simplex_with_noise');  
    [~, parms] = take_from_struct(parms, 'corr_type', 'Spearman');  %'Pearson',  'Spearman';
    [~, parms] = take_from_struct(parms, 'record_scores', true);  
    [~, parms] = take_from_struct(parms, 'rand_seed', 42);
    [~, parms] = take_from_struct(parms, 'num_restarts', 50);
    [~, parms] = take_from_struct(parms, 'H_lambda', 0.1);  
    [~, parms] = take_from_struct(parms, 'W_lambda', 0);  
    [~, parms] = take_from_struct(parms, 'log_transform', false);  
    [~, parms] = take_from_struct(parms, 'subsample_repeats', 5);
    [~, parms] = take_from_struct(parms, 'nmf_method', 'alsActiveSet');  
    [~, parms] = take_from_struct(parms, 'draw_log_scale', false);
    [~, parms] = take_from_struct(parms, 'init_type', 'random'); %'svd'
    [~, parms] = take_from_struct(parms, 'init_subtype', 'noise10');
    [~, parms] = take_from_struct(parms, 'random_init_spread', 1);
    [~, parms] = take_from_struct(parms, 'structure_filter', nan);
    [~, parms] = take_from_struct(parms, 'structre_type', 'relations_parent_level');
    [~, parms] = take_from_struct(parms, 'gene_subset', 'all');    
    [~, parms] = take_from_struct(parms, 'gene_okaty_filter', 'all');    

    [~, parms] = take_from_struct(parms, 'do_sep_init', true);
    [~, parms] = take_from_struct(parms, 'sep_H_lambda', 0);    
    
    
end
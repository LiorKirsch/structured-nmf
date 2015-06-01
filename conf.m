function parms = conf(parms)


parms.num_types = take_from_struct(parms, 'num_types', 3);  
parms.maxiter = take_from_struct(parms, 'maxiter', 500);  
parms.loglevel = take_from_struct(parms, 'loglevel', 0);  
parms.W_constraints = take_from_struct(parms, 'W_constraints', 'on_simplex_with_noise');  
  % 'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise';
parms.corr_type = take_from_struct(parms, 'corr_type', 'Spearman');  %'Pearson',  'Spearman';


parms.record_scores = take_from_struct(parms, 'record_scores', true);  
parms.rand_seed = take_from_struct(parms, 'rand_seed', 42);
% The answer to life the universe and everything
parms.num_restarts = take_from_struct(parms, 'num_restarts', 50);


parms.H_lambda = take_from_struct(parms, 'H_lambda', 0.1);  
parms.W_lambda = take_from_struct(parms, 'W_lambda', 0);  

parms.log_transform = false;
[subsample_repeats, parms] = take_from_struct(parms, 'subsample_repeats', 5);
parms.nmf_method = take_from_struct(parms, 'nmf_method', 'alsActiveSet');  
% alg_list = {'alsPinv', 'alsActiveSet', 'mm', 'alsBlockpivot','cjlin', 'prob'}; 

parms.draw_log_scale = false;
parms.do_sep_init = take_from_struct(parms, 'do_sep_init', true);  
parms.init_type = take_from_struct(parms, 'init_type', 'random');  %'svd' , 'random'
parms.structure_filter = take_from_struct(parms, 'structure_filter', nan);

end
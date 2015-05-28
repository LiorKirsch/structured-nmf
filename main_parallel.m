init;
  

parms = conf(parms);
%TODO( ? need this? show_results = false;



parms.gene_subset = 'all' ;%'okaty_infogain5000'; 'okaty_anova5000';% 'all';
parms.gene_okaty_filter = 'all'; % 'all'; 'cortex'; 'cortex_doyle'
% [gene_info, expression, parms] = gene_subset_selection(gene_info, ...
%                                                   expression, parms);


parms.structre_type = 'relations_parent_level';
parms = get_structure_matrix(parms.dataset_file, parms.structre_type,region_names, parms);
%=======================================================================


% % ====== load cell type specific genetic markers =======
% parms =load_markers(parms.dataset_file, gene_info, size(X{1},2), length(X), parms);

parms.do_sep_init = true;
parms.num_types = 3;
parms.W_constraints = 'on_simplex_with_noise';   % 'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise';
parms.nmf_method = 'alsActiveSet';
parms.num_restarts = 5; % <===  increase to 30
parms.W_lambda = 0;

parms.init_type = 'random';
% parms.num_markers = 20;

num_type_list = take_from_struct(parms, 'num_type_list', [3]) ;%1:8;
num_markers_list = take_from_struct(parms, 'num_markers_list', [5 20 50 100 200]);
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0 0.001 0.01 0.1 1 10 100 1000 inf]);

gene_subset_list = {'all','okaty_anova10000','okaty_anova5000'...
    'okaty_infogain5000','okaty_infogain10000'}; %,'okaty_anova1000'};
gene_okaty_filter_list = {'cortex_doyle','cortex','doyle', 'all'};

constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};


loop_over_var_name = {};
loop_over_var_value = {};
% loop_over_var_name{end + 1} = 'W_constraints';           % this cannot be the last list
% loop_over_var_value{end + 1} = constraints_list;       % this cannot be the last list
loop_over_var_name{end + 1} = 'gene_subset';           % this cannot be the last list
loop_over_var_value{end + 1} = gene_subset_list;       % this cannot be the last list
% loop_over_var_name{end + 1} = 'gene_okaty_filter';     % this cannot be the last list
% loop_over_var_value{end + 1} = gene_okaty_filter_list; % this cannot be the last list
% loop_over_var_name{end + 1} = 'num_types';
% loop_over_var_value{end + 1} = num_type_list;
% loop_over_var_name{end + 1} = 'num_markers';
% loop_over_var_value{end + 1} = num_markers_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;

parms.regions = region_names;
parms.cell_types = arrayfun(@(x) sprintf('#%d',x),1:parms.num_types,'UniformOutput',false);

results = loopOverHyperparms_real(X, gene_info, parms, loop_over_var_name, loop_over_var_value ,'');

parms.draw_log_scale = true;
draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
report_results(results);



% parms.H_lambda = 0.01;
% [cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
%     nmf(X, parms.num_types, 'alsWithRelations', parms);

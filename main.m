init;
parms = conf(parms);

dataset = take_from_struct(parms, 'dataset', 'brainspan2014');
[default_regions, parms.species] = get_region_set(dataset);
regions_short = take_from_struct(parms, 'regions_short', default_regions);
regions_short = sort(regions_short);

% Load all data
[parms.dataset_file expression, gross_region_vec, gene_info, ...
 gross_structures_info] = load_all(dataset, regions_short);

% limit to a set of specific genes
parms.gene_subset = 'all' ;% 'all'; 'okaty_infogain5000'; 'okaty_anova5000'; 'okaty_gainratio5000'
parms.gene_okaty_filter = 'all'; % 'all'; 'cortex'; 'doyle;'cortex_doyle'
% [gene_info, expression, parms] = gene_subset_selection(gene_info, ...
%                                                   expression, parms);

expression = change_to_linear_scale(expression);

gross_regions = gross_structures_info(gross_region_vec);
[X, region_names] = split_to_cell(expression, gross_regions);
clear('expression', 'gross_region_vec', 'gross_regions', ...
      'gross_structures_info');

fprintf('main: Get the structure matrix\n');
parms.regions = region_names;
parms.structre_type = 'relations_parent_level';
parms = get_structure_matrix(parms.dataset_file, parms.structre_type, ...
                                           region_names, parms);
%=======================================================================


% % ====== load cell type specific genetic markers =======
% parms =load_markers(parms.dataset_file, gene_info, size(X{1},2), length(X), parms);


% parms.num_markers = 20;

num_type_list = take_from_struct(parms, 'num_type_list', [3]) ;%1:8;
num_markers_list = take_from_struct(parms, 'num_markers_list', [5 20 50 100 200]);
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0, 10.^[-3:0.1:3], inf] );
gene_subset_list = take_from_struct(parms, 'gene_subset_list', ...
                                    {'all','okaty_anova10000', ...
                    'okaty_anova5000' 'okaty_infogain5000', ...
                    'okaty_infogain10000', 'okaty_infogain1000'});
gene_okaty_filter_list = take_from_struct(parms, 'gene_okaty_filter_list', ...
                                       {'cortex_doyle', 'cortex', 'doyle', 'all'});
constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};


loop_over_var_name = {};
loop_over_var_value = {};
% loop_over_var_name{end + 1} = 'W_constraints';           % this cannot be the last list
% loop_over_var_value{end + 1} = constraints_list;       % this cannot be the last list
% loop_over_var_name{end + 1} = 'gene_subset';           % this cannot be the last list
% loop_over_var_value{end + 1} = gene_subset_list;       % this cannot be the last list
% loop_over_var_name{end + 1} = 'gene_okaty_filter';     % this cannot be the last list
% loop_over_var_value{end + 1} = gene_okaty_filter_list; % this cannot be the last list
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;
% loop_over_var_name{end + 1} = 'num_markers';
% loop_over_var_value{end + 1} = num_markers_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;

parms.regions = region_names;
parms.cell_types = arrayfun(@(idx) sprintf('#%d',idx),1:parms.num_types,'UniformOutput',false);

fprintf('main: loop over hyper parameters\n');
results = loopOverHyperparms_real(X, gene_info, parms, loop_over_var_name, loop_over_var_value ,'');

if take_from_struct(parms, 'do_plot', true);
    parms.draw_log_scale = true;
    draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
    draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
    report_results(results);
end




% parms.H_lambda = 0.01;
% [cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
%     nmf(X, parms.num_types, 'alsWithRelations', parms);

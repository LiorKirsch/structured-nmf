init;
  
parms = conf(parms);


dataset = take_from_struct(parms, 'dataset', 'brainspan2014');
[default_regions, parms.species] = get_region_set(dataset);

regions_short = take_from_struct(parms, 'regions_short', default_regions);
regions_short = sort(regions_short);

switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions_short,'_'));
      [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
          = load_expression_and_regions('brainspan_rnaseq', regions_short);
    
    case 'zapala2005',    %==== Zapala selected regions ===
      parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('zapalaMouse', regions_short);
      gross_structures_info{strcmp(gross_structures_info, 'Isocortex')} ...
          = 'Cerebral_cortex';

      case 'human6' , %==== Human6 selected regions ===
        parms.dataset_file = 'Human6_selected_regions';
        [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
         ~] = load_expression_and_regions('human6LimitRegions', regions_short);
        gross_structures_info = gross_structures_info(:,4);
        gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x), ...
                                        gene_info.entrez_ids, ...
                                        'UniformOutput',false);
    otherwise 
      error('invalid dataset = [%s]\n', dataset);
end

% limit to a set of specific genes

parms.gene_subset = 'all' ;% 'all'; 'okaty_infogain5000'; 'okaty_anova5000'; 'okaty_gainratio5000'
parms.gene_okaty_filter = 'all'; % 'all'; 'cortex'; 'doyle;'cortex_doyle'
% [gene_info, expression, parms] = gene_subset_selection(gene_info, ...
%                                                   expression, parms);

expression = change_to_linear_scale(expression);

gross_regions = gross_structures_info(gross_region_vec);
[X,region_names] = split_to_cell(expression, gross_regions);
clear('expression', 'gross_region_vec', 'gross_regions', ...
      'gross_structures_info');

parms.structre_type = 'relations_parent_level';
parms = get_structure_matrix(parms.dataset_file, parms.structre_type,region_names, parms);
%=======================================================================


% % ====== load cell type specific genetic markers =======
% parms =load_markers(parms.dataset_file, gene_info, size(X{1},2), length(X), parms);


parms.W_constraints = 'on_simplex_with_noise';   % 'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise';
parms.init_type = 'random';

num_type_list = take_from_struct(parms, 'num_type_list', [3]) ;%1:8;
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0 0.001 0.01 0.1 1 10 100]);

gene_subset_list = take_from_struct(parms, 'gene_subset_list', {'all'});
% ,'okaty_anova10000', 'okaty_anova5000' 'okaty_infogain5000', 'okaty_infogain10000
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

results = loopOverResults_real(X, gene_info, parms, loop_over_var_name, loop_over_var_value ,'');

do_plot = take_from_struct(parms, 'do_plot', true);
if do_plot
    parms.draw_log_scale = true;
    draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
    draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
end
[best_score, base_score] = report_results(results, parms);

% parms.H_lambda = 0.01;
% [cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
%     nmf(X, parms.num_types, 'alsWithRelations', parms);

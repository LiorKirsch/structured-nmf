init;
  

parms = conf(parms);
%TODO( ? need this? show_results = false;

dataset = take_from_struct(parms, 'dataset', 'brainspan2014');
[default_regions, parms.species] = get_region_set(dataset);
regions = take_from_struct(parms, 'regions', default_regions);
regions = sort(regions);


switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions,'_'));
      [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
          = load_expression_and_regions('brainspan_rnaseq', regions);
    
    case 'zapala2005',    %==== Zapala selected regions ===
      parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('zapalaMouse', regions);
      gross_structures_info{strcmp(gross_structures_info, 'Isocortex')} ...
          = 'Cerebral_cortex';

      case 'human6' , %==== Human6 selected regions ===
        parms.dataset_file = 'Human6_selected_regions';
        [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
         ~] = load_expression_and_regions('human6LimitRegions', regions);
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

parms.regions = region_names;
parms.structre_type = 'relations_parent_level';
parms = get_structure_matrix(parms.dataset_file, parms.structre_type,region_names, parms);
%=======================================================================


% % ====== load cell type specific genetic markers =======
% parms =load_markers(parms.dataset_file, gene_info, size(X{1},2), length(X), parms);


parms.W_constraints = 'on_simplex_with_noise';   % 'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise';
parms.init_type = 'random'; %'svd' , 'random'
% parms.num_markers = 20;

num_type_list = take_from_struct(parms, 'num_type_list', [3]) ;%1:8;
num_markers_list = take_from_struct(parms, 'num_markers_list', [5 20 50 100 200]);
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0 0.001 0.01 0.1 1 10 100 1000 inf]);
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

results = loopOverHyperparms_real(X, gene_info, parms, loop_over_var_name, loop_over_var_value ,'');

parms.draw_log_scale = true;
draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
draw_indv_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
% report_results(results);



% parms.H_lambda = 0.01;
% [cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
%     nmf(X, parms.num_types, 'alsWithRelations', parms);

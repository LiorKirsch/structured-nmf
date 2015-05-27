init;
  

parms = conf(parms);

show_results = false;

dataset = take_from_struct(parms, 'dataset', 'brainspan2014');

switch dataset
    case 'kang2011',      %==== Kang ===
    
      parms.dataset_file = 'kang_regions';

      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
  case 'brainspan2014', %==== Brainspan ===
      [default_region_set, parms.species] = get_region_set('', dataset);
      subset_regions = take_from_struct(parms, 'subset_regions', default_region_set);
      
      parms.dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(subset_regions,'_'));
      [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
          = load_expression_and_regions('brainspan_rnaseq', subset_regions);
    
    case 'zapala2005',    %==== Zapala selected regions ===
      parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
      parms.species = 'mouse';
      regions_to_keep = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
                         'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
                         'Striatum';'Cerebellum';'Medulla'};
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('zapalaMouse', regions_to_keep);
      gross_structures_info{strcmp(gross_structures_info, 'Isocortex')} ...
          = 'Cerebral_cortex';

      case 'human6' , %==== Human6 selected regions ===
        parms.dataset_file = 'Human6_selected_regions';
        [regions_to_keep, parms.species] = get_region_set('', dataset);        
        [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
         ~] = load_expression_and_regions('human6LimitRegions', ...
                                          regions_to_keep);
        gross_structures_info = gross_structures_info(:,4);
        gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x), ...
                                        gene_info.entrez_ids, ...
                                        'UniformOutput',false);
    otherwise 
      error('invalid dataset = [%s]\n', dataset);
end

% limit to a set of specific genes

parms.gene_subset = 'all' ;%'okaty_infogain5000'; % 'all';
parms.gene_okaty_filter = 'all'; % 'all'; 'cortex'; 'cortex_doyle'
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

parms.do_sep_init = true;
parms.num_types = 3;
parms.W_constraints = 'on_simplex_with_noise';
parms.nmf_method = 'alsActiveSet';
parms.num_restarts = 10; % <===  increase to 30
parms.W_lambda = 0;

num_type_list = take_from_struct(parms, 'num_type_list', 3) ;%1:8;
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0 0.001 0.01 0.1 1 10 100 1000 inf]);

gene_subset_list = {'okaty_anova1000', 'all'};
gene_okaty_filter_list = {'cortex_doyle'};

loop_over_var_name = {};
loop_over_var_value = {};
loop_over_var_name{end + 1} = 'gene_subset';           % this cannot be the last list
loop_over_var_value{end + 1} = gene_subset_list;       % this cannot be the last list
loop_over_var_name{end + 1} = 'gene_okaty_filter';     % this cannot be the last list
loop_over_var_value{end + 1} = gene_okaty_filter_list; % this cannot be the last list
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;

parms.regions = region_names;
parms.cell_types = arrayfun(@(x) sprintf('#%d',x),1:parms.num_types,'UniformOutput',false);

results = loopOverHyperparms_real(X, gene_info, parms, loop_over_var_name, loop_over_var_value ,'');
report_results(results{1});



% parms.H_lambda = 0.01;
% [cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
%     nmf(X, parms.num_types, 'alsWithRelations', parms);

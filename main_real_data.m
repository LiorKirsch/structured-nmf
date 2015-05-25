init;
  

parms = conf(parms);

show_results = false;

dataset = take_from_struct(parms, 'dataset', 'brainspan2014');


% LIOR I SUGGEST THAT WE SET DATASETFILE TO BE STANDARD:
% parms.dataset_file = sprintf('%s_%s', dataset, regionset_name);
% The current setup is too manual and limits how many experiments
% we can run.  

%%% MY suggestion
% 1. Each region is assigned a single letter. eg 'M' for 'Motor
%    cortex'.Filenames and parameters operate based on this single letter 
% 2. Each region is also assigned a full name, this will be used in
%    printings and figures if necessary. 
% 3. region sets are either a set of letters, like MSV (would be
%    motor, sensory, visual cortices). or predefined names, like
%    'all'. or other sets. 
% 4. We define a function that returns the set of regions based on
%    regionset_name.   
%    region_set = get_region_set(regionset_name, dataset, species);
%    region_set = get_region_set(r'MSV', 'kang2011', 'human');




switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      parms.species = 'human';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.species = 'human';
      default_subset_regions =  {'A1C','AMY','CBC','DFC','HIP', ...
                          'IPC','ITC', 'M1C','MFC','OFC','S1C', ...
                          'STC','V1C','VFC'};   
      default_subset_regions =  {'DFC','OFC'}; 
      subset_regions = take_from_struct(parms, 'subset_regions', default_subset_regions);      
      
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
        parms.species = 'human';
        regions_to_keep = {'Frontal Lobe';'Cingulate gyrus';'hippocampal formation';...
                           'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';'Amygdala';'Basal Forebrain';...
                           'Striatum';'Thalamus';'Mesencephalon';'Cerebellar Cortex';...
                           'Myelencephalon'};
        % regions_to_keep = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
        %    'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
        %    'Striatum';'Cerebellum';'Medulla'};
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

parms.gene_subset = 'okaty_infogain5000' ;%'okaty_infogain5000'; % 'all';
parms.gene_okaty_filter = 'cortex_doyle'; % 'all';
[gene_info, expression, parms] = gene_subset_selection(gene_info, ...
                                                  expression, parms);

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
parms.num_restarts = 5; % <===  increase to 30
parms.W_lambda = 0;


num_type_list = 3 ;%1:8;
H_lambda_list = [0 0.001 0.01 0.1 1 10 100 1000 inf];

loop_over_var_name = {};
loop_over_var_value = {};
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

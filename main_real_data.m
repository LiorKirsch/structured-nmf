init;


parms = conf(parms);

show_results = false;

% % %==== Human6 selected regions ===
% % parms.dataset_file = 'Human6_selected_regions';
% % parms.species = 'human';
% % regions_to_keep = {'Frontal Lobe';'Cingulate gyrus';'hippocampal formation';'parahippocampal gyrus';...
% %                    'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';'Amygdala';'Basal Forebrain';...
% %                    'Striatum';'Hypothalamus';'Dorsal Thalamus';'Mesencephalon';'Cerebellar Cortex';...
% %                    'Pontine Tegmentum';'Myelencephalon'};
% % % regions_to_keep = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
% % %     'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
% % %     'Striatum';'Cerebellum';'Medulla'};
% % [expression, gross_region_vec, gene_info, ~, gross_structures_info, ~] = load_expression_and_regions('human6LimitRegions', regions_to_keep);
% gross_structures_info{strcmp(gross_structures_info,'Isocortex')} = 'Cerebral_cortex';

%==== Kang ===
% parms.dataset_file = 'kang_all_regions';
% parms.species = 'human';
% [expression, gross_region_vec, gene_info, ~, gross_structures_info, ~] = load_expression_and_regions('kangAllRegions', []);


% %==== Brainspan ===
% parms.dataset_file = 'brainspan_rnaseq';
% parms.species = 'human';
% [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
%     = load_expression_and_regions('brainspan_rnaseq', []);


% %==== Zapala selected regions ===
parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
parms.species = 'mouse';
regions_to_keep = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
    'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
    'Striatum';'Cerebellum';'Medulla'};
[expression, gross_region_vec, gene_info, ~, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', regions_to_keep);
gross_structures_info{strcmp(gross_structures_info,'Isocortex')} = 'Cerebral_cortex';


expression = change_to_linear_scale(expression);
gross_regions = gross_structures_info(gross_region_vec);
[X,region_names] = split_to_cell(expression, gross_regions);
clear('expression','gross_region_vec','gross_regions','gross_structures_info');


parms.structre_type = 'relations_parent_level';
parms = get_structure_matrix(parms.dataset_file, parms.structre_type,region_names, parms);
%=======================================================================


% % ====== load cell type specific genetic markers =======
% parms =load_markers(parms.dataset_file, gene_info, size(X{1},2), length(X), parms);


parms.do_sep_init = false;

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

parms.regions = mix_data.region;
parms.cell_types = mix_data.cell_types;

loopOverHyperparms_real(X, parms, loop_over_var_name, loop_over_var_value ,'');
                                 


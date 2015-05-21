%
% main_brainspan
%

init;

parms = conf(parms);

%==== Brainspan ===
[expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
    = load_expression_and_regions('brainspan_rnaseq', []);

gross_regions = gross_structures_info(gross_region_vec);

[all_X, all_region_names] = split_to_cell(expression, gross_regions);
clear('expression', 'gross_region_vec', 'gross_regions', ...
      'gross_structures_info');



% 

parms.relation_regions = {'DFC' 'OFC'}; % Regions to test
[~, relevant_regions] = ismember(parms.relation_regions, all_region_names);
X = all_X(relevant_regions);
region_names = all_region_names(relevant_regions);
parms.structure_matrix = [0 1; 1 0];

% parms.structre_type = 'relations_parent_level';
% [parms.structure_matrix,parms.relation_regions] = zapala_tree_structure(false,region_names);
% [parms.structure_matrix,parms.relation_regions] = ...
%    get_relation_structure(parms.structure_matrix, ...
%    parms.relation_regions, region_names, parms.structre_type);

parms.do_sep_init = false;
parms.num_types = 3;
parms.W_constraints = 'on_simplex_with_noise';
parms.nmf_method = 'alsActiveSet';
parms.num_restarts = 5; % <===  increase to 30
parms.W_lambda = 0;

fprintf('======== nmf seperate ========\n');
parms.H_lambda = 0;
% [cell_mix_sep.proportions, cell_mix_sep.cell_types] = nmf(X, parms.num_types, 'alsWithRelations', parms);  
[cell_mix_sep.proportions, cell_mix_sep.celltype_profile] = cellfun(@(x) ...
                                                  nmf(x,parms.num_types, ...
                                                  parms.nmf_method, ...
                                                  parms), X, ...
                                                  'UniformOutput', false);

fprintf('======== nmf single ========\n');
parms.H_lambda = inf;
[cell_mix_single.proportions, cell_mix_single.celltype_profile] = ...
    nmf(X, parms.num_types, 'alsWithRelations', parms);
cell_mix_single2.cell_types = {'type #1', 'type #2', 'type #3'};
cell_mix_single2.celltype_profile = cell_mix_single.celltype_profile{1}';

[cell_mix_sep2.celltype_profile, cell_mix_sep2.cell_types] = ...
    join_profiles(cell_mix_sep.celltype_profile, region_names);

baseline.celltype_profile = cellfun(@(x) repmat(mean(x),3,1), ...
                                    X,'UniformOutput',false);
baseline.proportions = cellfun(@(x) zeros(3,3), X ,'UniformOutput', false);

% [baseline.celltype_profile, ~] =
% join_profiles(baseline.celltype_profile, region_names);
% baseline.celltype_profile = baseline.celltype_profile';
% baseline.cell_types = region_names;

compare_to_true_profile( cellfun(@(x) x', ...
                                 cell_mix_sep.celltype_profile, ...
                                 'UniformOutput',false), ...
                         cell_mix_sep.proportions, gene_info, ...
                         region_names,'human',parms);
compare_to_true_profile( cellfun(@(x) x', ...
                                 cell_mix_single.celltype_profile, ...
                                 'UniformOutput',false), ...
                         cell_mix_single.proportions, gene_info, ...
                         region_names,'human',parms);
compare_to_true_profile( cellfun(@(x) x', ...
                                 baseline.celltype_profile, ...
                                 'UniformOutput',false), ...
                         baseline.proportions, gene_info, ...
                         region_names,'human',parms);

% compare celltype profile with celltype profiles from Okaty PLoS One
parms.species = 'human';
cell_mix_sep2.celltype_profile = cell_mix_sep2.celltype_profile';
[aucs,median_within, median_outside] = compare_nmf_to_doyle(cell_mix_sep2, gene_info, parms);
[aucs,median_within, median_outside] = compare_nmf_to_doyle(cell_mix_single2, gene_info, parms);
[aucs,median_within, median_outside] = compare_nmf_to_doyle(baseline, gene_info, parms);
title(dataset_name);

draw_proprtions_for_regions(cell_mix.proportions, cell_mix.cell_types, ...
                            gross_structures_info, gross_region_vec);
ylim([0,0.8]);
show_proportions(cell_mix.proportions, cell_mix.cell_types);

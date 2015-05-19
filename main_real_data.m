init;
addpath('/home/lab/lior/Projects/load datasets/');
addpath('wet celltype compare/');

parms = conf(parms);
%===  Zapala data  ===
% dataset_name = 'zapala';
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', []);
% gross_structures_info(strcmp('Bed nucleus of the stria terminalis',gross_structures_info))  = {'BNST'};


%==== Zapala cortex and hippo ===
% dataset_name = 'zapala cortex';
% regions_to_keep = {'Olfactory Bulbs';'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';};
% [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', regions_to_keep);


%==== Zapala cortex and hippo ===
dataset_name = 'Zapala isocortex, medulla, striatum and cerebellum';
regions_to_keep = {'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
    'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
    'Striatum';'Cerebellum';'Medulla'};
[expression, gross_region_vec, gene_info, ~, gross_structures_info, ~] = load_expression_and_regions('zapalaMouse', regions_to_keep);
gross_structures_info{strcmp(gross_structures_info,'Isocortex')} = 'Cerebral_cortex';


% mouse_cell_types = load('mouse_cell_type_profiles.mat');
% [~, reorder_predicted] = reorderUsingId(mouse_cell_types.all_symbols, gene_info.gene_symbols);
% expression = expression(:,reorder_predicted);
% gene_info.gene_symbols = gene_info.gene_symbols(reorder_predicted);
% gene_info.entrez_ids = gene_info.entrez_ids(reorder_predicted);
% gene_info.gene_full_name = gene_info.gene_full_name(reorder_predicted);
% gene_info.probe_id = gene_info.probe_id(reorder_predicted);


gross_regions = gross_structures_info(gross_region_vec);
[X,region_names] = split_to_cell(expression, gross_regions);
clear('expression','gross_region_vec','gross_regions','gross_structures_info');
parms.structre_type = 'relations_parent_level';
[parms.structure_matrix,parms.relation_regions] = zapala_tree_structure(false,region_names);
[parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,region_names,parms.structre_type);
%     figure;imagesc(parms.structure_matrix);colorbar; colormap(jet);
%     ax = gca;
%     ax.XTick = 1:length(parms.relation_regions);
%     ax.XTickLabel = parms.relation_regions;
%     ax.YTick = 1:length(parms.relation_regions);
%     ax.YTickLabel = parms.relation_regions;
%     ax.XTickLabelRotation	=45;

%=======================================================================


parms.do_sep_init = false;

parms.num_types = 3;
parms.W_constraints = 'on_simplex_with_noise';
parms.nmf_method = 'alsActiveSet';
parms.num_restarts = 30; % <===  increase to 30
parms.W_lambda = 0;


neuro_markers = {'Stmn2','Celf4','Syt1','Gria3','Dlg3','Dlg4','Tubb3','Map2'};
astro_markers = { 'Gfap' , 'Aqp4' ,'Fgfr3' ,'Slc1a2' ,'Gjb6' };
oligo_markers = {'Mbp' ,'Sox10' ,'Mag' ,'Mog'};
H_markers = false(parms.num_types,size(X{1},2));
H_markers(1, ismember(gene_info.gene_symbols,neuro_markers) ) = true;
H_markers(2, ismember(gene_info.gene_symbols,astro_markers) ) = true;
H_markers(3, ismember(gene_info.gene_symbols,oligo_markers) ) = true;
parms.H_markers = H_markers;


fprintf('======== nmf seperate ========\n');
parms.H_lambda = 0;
% [cell_mix_sep.proportions, cell_mix_sep.cell_types] = nmf(X, parms.num_types, 'alsWithRelations', parms);  
[cell_mix_sep.proportions, cell_mix_sep.celltype_profile]=cellfun(@(x) ...
       nmf(x,parms.num_types, parms.nmf_method, parms) ,X,'UniformOutput',false);

fprintf('======== nmf single ========\n');
parms.H_lambda = inf;
[cell_mix_single.proportions, cell_mix_single.celltype_profile] = nmf(X, parms.num_types, 'alsWithRelations', parms);  
cell_mix_single2.cell_types = {'type #1', 'type #2', 'type #3'};
cell_mix_single2.celltype_profile = cell_mix_single.celltype_profile{1}';

[cell_mix_sep2.celltype_profile, cell_mix_sep2.cell_types] = join_profiles(cell_mix_sep.celltype_profile, region_names);

baseline.celltype_profile = cellfun(@(x) repmat(mean(x),3,1), X,'UniformOutput',false);
baseline.proportions = cellfun(@(x) zeros(3,3), X ,'UniformOutput',false);
% [baseline.celltype_profile, ~] = join_profiles(baseline.celltype_profile, region_names);
% baseline.celltype_profile = baseline.celltype_profile';
% baseline.cell_types = region_names;


compare_to_true_profile( cellfun(@(x) x',cell_mix_sep.celltype_profile,'UniformOutput',false),...
    cell_mix_sep.proportions, gene_info,region_names,parms);
compare_to_true_profile( cellfun(@(x) x',cell_mix_single.celltype_profile,'UniformOutput',false),...
    cell_mix_single.proportions, gene_info,region_names,parms);
compare_to_true_profile( cellfun(@(x) x',baseline.celltype_profile,'UniformOutput',false),...
    baseline.proportions, gene_info,region_names,parms);

% compare celltype profile with celltype profiles from Okaty PLoS One
parms.species = 'mouse';
cell_mix_sep2.celltype_profile = cell_mix_sep2.celltype_profile';
[aucs,median_within, median_outside] = compare_nmf_to_doyle(cell_mix_sep2, gene_info, parms);
[aucs,median_within, median_outside] = compare_nmf_to_doyle(cell_mix_single2, gene_info, parms);
[aucs,median_within, median_outside] = compare_nmf_to_doyle(baseline, gene_info, parms);
title(dataset_name);

draw_proprtions_for_regions(cell_mix.proportions, cell_mix.cell_types, gross_structures_info, gross_region_vec); ylim([0,0.8]);
show_proportions(cell_mix.proportions, cell_mix.cell_types);

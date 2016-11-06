parms.dummy = 1;
parms = conf(parms);

dataset = take_from_struct(parms, 'dataset', 'kang2011');

switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      parms.species = 'human';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.dataset_file = 'brainspan_rnaseq_DFC_OFC';
      parms.species = 'human';
      subset_regions = {'DFC' 'OFC'}; % Regions to test
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
      gross_structures_info{strcmp(gross_structures_info,'Isocortex')} = 'Cerebral_cortex';

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
        gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x),gene_info.entrez_ids,'UniformOutput',false);


    otherwise 
      error('invalid dataset = [%s]\n', dataset);
end
      
% limit to a set of specific genes
parms.gene_subset = 'all';
[gene_info, expression, parms] = gene_subset_selection(gene_info, expression, parms);

gross_regions = gross_structures_info(gross_region_vec);
[X,region_names] = split_to_cell(expression, gross_regions);


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
parms.H_lambda = 0.01;

[cell_mix.proportions, cell_mix.celltype_profile] = load_nmf_results(...
     X, parms.num_types, 'alsWithRelations', parms);

[cell_mix.celltype_profile, cell_mix.cell_types] = join_profiles(cell_mix.celltype_profile, region_names); 
cell_mix.celltype_profile = cell_mix.celltype_profile';
% cell_mix.celltype_profile = cellfun(@transpose, cell_mix.celltype_profile,'UniformOutput',false);
figure('Name','Seperate');compare_nmf_to_doyle(cell_mix, gene_info, parms);



baseline.celltype_profile = cellfun(@(x) mean(x), X,'UniformOutput',false);
baseline.proportions = cellfun(@(x) zeros(3,3), X ,'UniformOutput',false);
[baseline.celltype_profile, ~] = join_profiles(baseline.celltype_profile, region_names);
baseline.celltype_profile = baseline.celltype_profile';
baseline.cell_types = region_names;

figure('Name','Mean profile');compare_nmf_to_doyle(baseline, gene_info, parms);

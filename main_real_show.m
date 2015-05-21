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

baseline.celltype_profile = cellfun(@(x) mean(x), X,'UniformOutput',false);
baseline.proportions = cellfun(@(x) zeros(3,3), X ,'UniformOutput',false);
[baseline.celltype_profile, ~] = join_profiles(baseline.celltype_profile, region_names);
baseline.celltype_profile = baseline.celltype_profile';
baseline.cell_types = region_names;

% compare celltype profile with celltype profiles from Doyle

cell_mix_sep2.celltype_profile = cell_mix_sep2.celltype_profile';
figure('Name','Seperate');compare_nmf_to_doyle(cell_mix_sep2, gene_info, parms);
figure('Name','Unified');compare_nmf_to_doyle(cell_mix_single2, gene_info, parms);
figure('Name','Mean profile');compare_nmf_to_doyle(baseline, gene_info, parms);



% compare_to_true_profile( cellfun(@(x) x',cell_mix_sep.celltype_profile,'UniformOutput',false),...
%     cell_mix_sep.proportions, gene_info,region_names,'human',parms);
% compare_to_true_profile( cellfun(@(x) x',cell_mix_single.celltype_profile,'UniformOutput',false),...
%     cell_mix_single.proportions, gene_info,region_names,'human',parms);
% compare_to_true_profile( cellfun(@(x) x',baseline.celltype_profile,'UniformOutput',false),...
%     baseline.proportions, gene_info,region_names,'human',parms);

draw_proprtions_for_regions(cell_mix.proportions, cell_mix.cell_types, gross_structures_info, gross_region_vec); ylim([0,0.8]);
show_proportions(cell_mix.proportions, cell_mix.cell_types);
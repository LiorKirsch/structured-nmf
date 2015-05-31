parms.a=1;
parms = conf(parms);


dataset = take_from_struct(parms, 'dataset', 'human6');
dataset_mouse = take_from_struct(parms, 'dataset_mouse', 'zapala');

[regions, parms.species] = get_region_set(dataset);
regions = sort(regions);

switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      [human_expression, human_gross_region_vec, human_gene_info, ~, human_gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions,'_'));
      regions = {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
                         'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};
      [human_expression, human_gross_region_vec, human_gene_info, ~, human_gross_structures_info] ...
          = load_expression_and_regions('brainspan_rnaseq', regions);
    
   

      case 'human6' , %==== Human6 selected regions ===
        parms.dataset_file = 'Human6_selected_regions';
        [human_expression, human_gross_region_vec, human_gene_info, ~, human_gross_structures_info, ...
         ~] = load_expression_and_regions('human6LimitRegions', regions);
        human_gross_structures_info = human_gross_structures_info(:,4);
        human_gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x), ...
                                        human_gene_info.entrez_ids, ...
                                        'UniformOutput',false);
    otherwise 
      error('invalid dataset = [%s]\n', dataset);
end


switch dataset_mouse
    case 'zapala',      %==== Zapala ===
        [regions, parms.species] = get_region_set('zapala2005');
        regions = sort(regions);
        parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
          [mouse_expression, mouse_gross_region_vec, mouse_gene_info, ~, mouse_gross_structures_info, ...
           ~] = load_expression_and_regions('zapalaMouse', regions);
          mouse_gross_structures_info{strcmp(mouse_gross_structures_info, 'Isocortex')} ...
              = 'Cerebral_cortex';
    case 'doyle'
    otherwise 
      error('invalid dataset = [%s]\n', dataset_mouse);
end
  
  
  % ====== SORT THE GENE WHICH CAN BE MAPPED HOMOLOGOUS 1-1
  
  addpath('/cortex/code/matlab/homologous_gene_mapping/');
           
           
  [gene_to_group_mouse, gene_to_group_primate, homologous_group_id] =  gene_to_homolog_group(...
      'mouse_laboratory','human', mouse_gene_info.gene_symbols, 'symbol'...
                            ,human_gene_info.entrez_ids,'entrez_gene_ID');

    groups_with_1_to_1 = sum(gene_to_group_mouse,1) == 1  & sum(gene_to_group_primate,1) == 1;
    gene_to_group_mouse = gene_to_group_mouse(:,groups_with_1_to_1);
    gene_to_group_primate = gene_to_group_primate(:,groups_with_1_to_1);
    gene_to_group_mouse = (1:size(gene_to_group_mouse,1)) * gene_to_group_mouse ;
    gene_to_group_primate = (1:size(gene_to_group_primate,1)) * gene_to_group_primate ;

    mouse_expression = mouse_expression(:,gene_to_group_mouse);
    human_expression = human_expression(:,gene_to_group_primate);
    
    num_genes = length(gene_to_group_primate);
    fprintf('Computing correlation using %d genes which can be mapped between the datasets\n',num_genes);

  % ====== for each region calc the mean profile
  num_regions_mouse = length(mouse_gross_structures_info);
  mean_exp_mouse = nan(num_regions_mouse,num_genes);
  for i =1:num_regions_mouse
      mean_exp_mouse(i,:) = mean(mouse_expression( i==mouse_gross_region_vec,:),1);
  end
  
  num_regions_human = length(human_gross_structures_info);
  mean_exp_human = nan(num_regions_human,num_genes);
  for i =1:num_regions_human
      mean_exp_human(i,:) = mean(human_expression( i==human_gross_region_vec,:),1);
  end
  
  % ===== For each pair of regions calculate the calibrated correlation.
  disp('computing the callibrated correlation between the datasets');
  corr_matrix = callibrated_corr(mean_exp_mouse', mean_exp_human', 'pearson');
  imagesc(corr_matrix);
  colorbar; colormap(jet);
    ax = gca;
    ax.YTick = 1:length(mouse_gross_structures_info);
    ax.YTickLabel = mouse_gross_structures_info;
    ax.XTick = 1:length(human_gross_structures_info);
    ax.XTickLabel = human_gross_structures_info;
    ax.XTickLabelRotation	=45;
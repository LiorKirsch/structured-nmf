function compare_to_true_profile(predicted_profiles, predicted_proportions, gene_info, region_name,parms)

    % load doyle
    mouse_cell_types = load_doyle_true_type();
    
    % select gene which appear in both datasets according to gene symbol
    [reorder_mouse_cell_type, reorder_predicted] = reorderUsingId(mouse_cell_types.all_symbols, gene_info.gene_symbols);
    reorder_mouse_cell_type = mouse_cell_types.refer_to_index(reorder_mouse_cell_type);
    mouse_cell_types.expression = mouse_cell_types.expression(reorder_mouse_cell_type,:);
    predicted_profiles = cellfun(@(x) x(reorder_predicted,:), predicted_profiles,'UniformOutput',false);
    
    % get true profile for each region
    true_profiles = match_region_with_true_profile(mouse_cell_types, region_name);

    % get the score
    score = get_mean_score(predicted_profiles, predicted_proportions, true_profiles, region_name,parms);
end

function best_score = get_mean_score(predicted_profiles, predicted_proportions, true_profiles, region_name,parms)
num_regions = length(predicted_profiles);
assert(length(true_profiles) == num_regions,'true and predicted profiles should have the same number of elements');

    
best_scores = nan(num_regions,1);
for i = 1:num_regions
    % need to transpose the expression
    GT_proportions = zeros(size(predicted_proportions{i},1),size(true_profiles{i},2));
   [~, ~, best_scores(i), ~] = match_profiles_to_gt(predicted_proportions{i}, predicted_profiles{i}',...
       true_profiles{i}', GT_proportions', 'spearman'); 
   fprintf('%s - %g\n',region_name{i}, best_scores(i));
end

 best_score = mean_corr_coeff(best_scores);
 fprintf('====MEAN SCORE: %g====\n', best_score);
end

function celltype_expression = match_region_with_true_profile(mouse_cell_types, regions)
    celltype_expression = cell(length(regions),1);
    for i = 1:length(regions)
        curr_region = regions{i};
        
%         {'Cerebral_cortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
%     'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3';...
%     'Striatum';'Cerebellum';'Medulla'};

%         {'H_lva_cortex';'H_lvb_cortex';'H_l6_cortex';'H_cck+_cortex';'H_mn_brainstem';'H_chol_basal_fb';...
%             'H_chol_spinal';'H_chol_corp_striatum';'H_cort+_cortex';'H_drd1+_msn_striatum';'H_drd2+_msn_striatum';...
%             'H_golgi_cerebellum';'H_uni_brush_cerebellum';'H_stell_basket_cerebellum';'H_granule_cerebellum';...
%             'H_molig_cerebellum';'H_molig_cortex';'H_mixed_olig_cerebellum';'H_mixed_olig_cortex';...
%             'H_purkinje_cerebellum';'H_pnoc+_cortex';'H_bergmann_cerebellum';'H_astro_cerebellum';'H_astro_cortex'}

        switch curr_region
            case {'Cerebral_cortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';...
                  'Hippocampus';'Dentate Gyrus';'Hippocampus CA1';'Hippocampus CA3'}
            
                neuro_id = {'H_lva_cortex';'H_lvb_cortex';'H_l6_cortex';'H_cck+_cortex';...
                            'H_cort+_cortex'; 'H_pnoc+_cortex'};
                oligo_id = {'H_molig_cortex';'H_mixed_olig_cortex'};
                astro_id = {'H_astro_cortex'};
            case 'Striatum'
                neuro_id = { 'H_chol_corp_striatum';'H_drd1+_msn_striatum';'H_drd2+_msn_striatum'} ;
                oligo_id = {};
                astro_id = {};
            case 'Cerebellum'
                neuro_id = { 'H_golgi_cerebellum';'H_stell_basket_cerebellum';'H_granule_cerebellum';...
                    'H_purkinje_cerebellum'};
                oligo_id = {'H_molig_cerebellum';'H_mixed_olig_cerebellum'};
                astro_id = {'H_uni_brush_cerebellum';'H_bergmann_cerebellum';'H_astro_cerebellum'};
            case 'Medulla'
                neuro_id = { 'H_mn_brainstem'};
                oligo_id = {};
                astro_id = {};
            otherwise
                error('unknown region %s' ,curr_region)
        end
        
        expression = [];
        filterneuron = ismember(mouse_cell_types.cell_type_id, neuro_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
  
        filterneuron = ismember(mouse_cell_types.cell_type_id, astro_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
        
        filterneuron = ismember(mouse_cell_types.cell_type_id, oligo_id);
        sample2type = mouse_cell_types.sample2type(:,filterneuron);
        if ~isempty(sample2type)
            all_samples = any(sample2type,2);
            mean_expression = mouse_cell_types.expression * double(all_samples) / sum(all_samples);
            expression = cat(2,expression, mean_expression);
        end
        
        celltype_expression{i} = expression;
    end


end

function [expressionA,expressionB] = limit_to_genes_in_intersection(expressionA,expressionB, gene_list_A, gene_list_B)

    [reorder_A, reorder_B] = reorderUsingId(gene_list_A, gene_list_B);
    reorder_mouse_cell_type = mouse_cell_types.refer_to_index(reorder_A);

    expressionA = expressionA(reorder_mouse_cell_type,:);
    expressionB = expressionB(reorder_B,:);
    
end

function mouse_cell_types = load_doyle_true_type()


    mouse_cell_types = load('mouse_cell_type_profiles.mat');
    mouse_cell_types.expression = 2.^mouse_cell_types.expression;


    cell_type_filter = strmatch('Doyle', mouse_cell_types.reference);
    mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, cell_type_filter);
    
    exp = mouse_cell_types.expression * double(mouse_cell_types.sample2type);
    mouse_cell_types.expression_types = exp ./repmat(sum(mouse_cell_types.sample2type,1),size(mouse_cell_types.expression,1),1);
  
end
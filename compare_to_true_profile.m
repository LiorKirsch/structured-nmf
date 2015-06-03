function [gene_inds_true_type,gene_inds_predictions, ...
          gene_symb_true_type, gene_symb_predictions] = ...
        compare_to_true_profile(mouse_cell_types,gene_info, species,parms)

    % load doyle
%     mouse_cell_types = load_doyle_true_type();
    
    switch species
        case 'mouse'
            
            % select gene which appear in both datasets according to gene symbol
            [reorder_mouse_cell_type, reorder_predicted] = ...
                reorderUsingId(mouse_cell_types.all_symbols, ...
                               gene_info.gene_symbols);
            reorder_mouse_cell_type = ...
                mouse_cell_types.refer_to_index(reorder_mouse_cell_type);
            % mouse_cell_types.expression =
            % mouse_cell_types.expression(reorder_mouse_cell_type,:);
            % predicted_profiles = cellfun(@(x)
            % x(reorder_predicted,:),
            % predicted_profiles,'UniformOutput',false);
            gene_inds_true_type = reorder_mouse_cell_type;
            gene_inds_predictions = reorder_predicted;
            
            gene_symb_true_type  = mouse_cell_types.gene_symbol(reorder_mouse_cell_type);
            gene_symb_predictions  = gene_info.gene_symbols(reorder_predicted);
            
      case {'monkey', 'human'}
            addpath('/cortex/code/matlab/homologous_gene_mapping/');
            switch species
                case 'monkey'
                    [gene_to_group_mouse, gene_to_group_primate, ...
                     homologous_group_id] = ...
                        gene_to_homolog_group('mouse_laboratory', ...
                                              'rhesus_macaque', ...
                                              mouse_cell_types.gene_symbol, ...
                                              'symbol', gene_info ...
                                              .entrez_ids, 'entrez_gene_ID');
                case 'human'
                    [gene_to_group_mouse, gene_to_group_primate, ...
                     homologous_group_id] = ...
                        gene_to_homolog_group('mouse_laboratory','human', ...
                                              mouse_cell_types.gene_symbol, ...
                                              'symbol', gene_info ...
                                              .entrez_ids, 'entrez_gene_ID');
            end
            
            groups_with_1_to_1 = sum(gene_to_group_mouse, 1) == 1  & sum(gene_to_group_primate,1) == 1;
            gene_to_group_mouse = gene_to_group_mouse(:,groups_with_1_to_1);
            gene_to_group_primate = gene_to_group_primate(:,groups_with_1_to_1);
            gene_to_group_mouse = (1:size(gene_to_group_mouse,1)) * gene_to_group_mouse ;
            gene_to_group_primate = (1:size(gene_to_group_primate,1)) * gene_to_group_primate ;

            gene_inds_true_type = gene_to_group_mouse;
            gene_inds_predictions = gene_to_group_primate;
            
            gene_symb_true_type  = mouse_cell_types.gene_symbol(gene_to_group_mouse);
            gene_symb_predictions  = gene_info.gene_symbols(gene_to_group_primate);
    end
    
    % get true profile for each region
%     true_profiles = match_region_with_true_profile(mouse_cell_types, region_name);

    % get the score
%     score = get_mean_score(predicted_profiles, predicted_proportions, true_profiles, region_name,parms);
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
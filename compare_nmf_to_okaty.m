function [aucs,median_within, median_outside] = compare_nmf_to_okaty(cell_mix, gene_info, parms)


    mouse_cell_types = load('mouse_cell_type_profiles.mat');
    % [found_translation, entrez_translation] = translateAffyToEntrez(mouse_cell_types.gene_affi_id, 'entrez', false,gene_info.entrez_ids);
    % mouse_cell_types.gene_entrez = entrez_translation;
    % mouse_cell_types.gene_symbol = mouse_cell_types.gene_symbol(found_translation);
    % 
    % [found_translation, symbol_translation] = translateAffyToEntrez(mouse_cell_types.gene_affi_id, 'symbol', false,gene_info.gene_symbols);
    % mouse_cell_types.gene_symbol = symbol_translation;
    % 
    % mouse_cell_types.expression = mouse_cell_types.expression(found_translation,:);
    % mouse_cell_types.gene_affi_id = mouse_cell_types.gene_affi_id(found_translation);
    
    if ~exist('limit_to_cortical_cell_types','var');
        parms.limit_to_cortical_cell_types = false;
    end
    [mouse_cell_types,seperation_indices, cell_type_order] = order_sample_by_type(mouse_cell_types, parms.limit_to_cortical_cell_types);
    
    switch parms.species
        case 'mouse'
            [reorder_mouse_cell_type, reorder_zapala] = reorderUsingId(mouse_cell_types.all_symbols, gene_info.gene_symbols);
            reorder_mouse_cell_type = mouse_cell_types.refer_to_index(reorder_mouse_cell_type);
            reorder_mouse_cell_type_expresion = mouse_cell_types.expression(reorder_mouse_cell_type,:);
            reorder_cellmix_expresion = cell_mix.celltype_profile(reorder_zapala,:);
        case {'monkey', 'human'}
            addpath('/cortex/code/matlab/homologous_gene_mapping/');
            switch parms.species
                case 'monkey'
                    [gene_to_group_mouse, gene_to_group_primate, homologous_group_id] =  gene_to_homolog_group('mouse_laboratory','rhesus_macaque', mouse_cell_types.gene_symbol, 'symbol',gene_info.entrez_ids,'entrez_gene_ID');
                case 'human'
                    [gene_to_group_mouse, gene_to_group_primate, homologous_group_id] =  gene_to_homolog_group('mouse_laboratory','human', mouse_cell_types.gene_symbol, 'symbol',gene_info.entrez_ids,'entrez_gene_ID');
            end
            
            groups_with_1_to_1 = sum(gene_to_group_mouse,1) == 1  & sum(gene_to_group_primate,1) == 1;
            gene_to_group_mouse = gene_to_group_mouse(:,groups_with_1_to_1);
            gene_to_group_primate = gene_to_group_primate(:,groups_with_1_to_1);
            gene_to_group_mouse = (1:size(gene_to_group_mouse,1)) * gene_to_group_mouse ;
            gene_to_group_primate = (1:size(gene_to_group_primate,1)) * gene_to_group_primate ;

            reorder_mouse_cell_type_expresion = mouse_cell_types.expression(gene_to_group_mouse,:);
            reorder_cellmix_expresion = cell_mix.celltype_profile(gene_to_group_primate,:);
    end

    fprintf('Computing correlation using %d genes which can be mapped between the datasets\n',size(reorder_mouse_cell_type_expresion,1) );
    corr_matrix = corr(reorder_mouse_cell_type_expresion, reorder_cellmix_expresion, 'type','spearman');
    draw_imagesc_with_seperation(corr_matrix,seperation_indices,cell_type_order);

    ax = gca;
    ax.XTick = 1:length(cell_mix.cell_types);
    ax.XTickLabel = cell_mix.cell_types;
    ax.XTickLabelRotation	=45;
    % ax.YTick = 1:length(mouse_cell_types.cellTypesDescription);
    % ax.YTickLabel = mouse_cell_types.cellTypesDescription;
    
    full_name = set_filenames('figure_confusion', parms);
    saveas(gca,[full_name,'_corr'],'png');
    figure;
    single_sample_scatter(reorder_mouse_cell_type_expresion(:,28), reorder_cellmix_expresion(:,strcmp('Neurons',cell_mix.cell_types)),parms);
    saveas(gca,[full_name,'_scatter'],'png');
    sample_cell_type_id =double(mouse_cell_types.sample2type) * ( ( 1:size(mouse_cell_types.sample2type,2))');
    mouse_cell_types.is_neuron = mouse_cell_types.is_neuron(sample_cell_type_id);
    mouse_cell_types.is_astro = mouse_cell_types.is_astro(sample_cell_type_id);
    mouse_cell_types.is_oligo = mouse_cell_types.is_oligo(sample_cell_type_id);
    aucs = drawROC(cell_mix, mouse_cell_types, corr_matrix);
    saveas(gca,[full_name,'_roc'],'png');
    [median_within, median_outside] = compute_median_corr(corr_matrix, mouse_cell_types);
    
end

function [median_within, median_outside] = compute_median_corr(corr_matrix, mouse_cell_types)
    astro_median_same = median( corr_matrix(mouse_cell_types.is_astro,1));
    astro_median_diff = median( corr_matrix(~mouse_cell_types.is_astro,1));
    
    neuron_median_same = median( corr_matrix(mouse_cell_types.is_neuron,2));
    neuron_median_diff = median( corr_matrix(~mouse_cell_types.is_neuron,2));
    
    oligo_median_same = median( corr_matrix(mouse_cell_types.is_oligo,3));
    oligo_median_diff = median( corr_matrix(~mouse_cell_types.is_oligo,3));
    
    median_within = [astro_median_same, neuron_median_same, oligo_median_same];
    median_outside = [astro_median_diff, neuron_median_diff, oligo_median_diff];
end


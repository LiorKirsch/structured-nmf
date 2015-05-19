function [aucs,median_within, median_outside] = compare_nmf_to_doyle(cell_mix, gene_info, parms)


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
    mouse_cell_types.expression = 2.^mouse_cell_types.expression;
    if ~exist('limit_to_cortical_cell_types','var');
        parms.limit_to_cortical_cell_types = false;
    end

    [mouse_cell_types,seperation_indices, cell_type_order] = order_sample_by_type(mouse_cell_types, parms.limit_to_cortical_cell_types);
    
    cell_type_filter = strmatch('Doyle', mouse_cell_types.reference);
    mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, cell_type_filter);
    exp = mouse_cell_types.expression * double(mouse_cell_types.sample2type);
    exp = exp ./repmat(sum(mouse_cell_types.sample2type,1),size(mouse_cell_types.expression,1),1);
    
    switch parms.species
        case 'mouse'
            % since some genes have more than one name we first use
            % all_symbols which hold the full symbol list than we map it
            % back to the indices that are of intrests to us.
            [reorder_mouse_cell_type, reorder_zapala] = reorderUsingId(mouse_cell_types.all_symbols, gene_info.gene_symbols);
            reorder_mouse_cell_type = mouse_cell_types.refer_to_index(reorder_mouse_cell_type);
%             reorder_mouse_cell_type_expresion = mouse_cell_types.expression(reorder_mouse_cell_type,:);
% tmp_name_A = mouse_cell_types.gene_symbol(reorder_mouse_cell_type); %CHECKING
% tmp_name_B =  gene_info.gene_symbols(reorder_zapala); %CHECKING
            reorder_mouse_cell_type_expresion = exp(reorder_mouse_cell_type,:);
            reorder_cellmix_expresion = cell_mix.celltype_profile(reorder_zapala,:);
        otherwise
            error('unkown specices');
    end

    fprintf('Computing correlation using %d genes which can be mapped between the datasets\n',size(reorder_mouse_cell_type_expresion,1) );
%     corr_matrix = corr(reorder_mouse_cell_type_expresion, reorder_cellmix_expresion, 'type','pearson');
    corr_matrix = callibrated_corr(reorder_mouse_cell_type_expresion, reorder_cellmix_expresion, 'spearman');
    imagesc(corr_matrix);  
    colormap(jet); 
    imagescwithnan(corr_matrix,jet,[0 1 1]) %# [0 1 1] is cyan

    ax = gca;
    ax.XTick = 1:length(cell_mix.cell_types);
    ax.XTickLabel = cell_mix.cell_types;
    ax.XTickLabelRotation	=45;
    ax.YTick = 1:length(mouse_cell_types.cell_type_description);
    type_and_region = cellfun(@(x,y) sprintf('%s (%s)',x,y), mouse_cell_types.cell_type_description, mouse_cell_types.anatomical_region, 'UniformOutput',false);
    ax.YTickLabel = type_and_region;
    
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


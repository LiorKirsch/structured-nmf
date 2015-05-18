function [aucs,median_within, median_outside] = compare_deconv_to_giger_cell_type(cell_mix, gene_info, species,limit_to_cortical_cell_types)


    addpath('markers profiles/');
    giger_cell_types = load('giger_cell_profiles.mat');
    
    if ~exist('limit_to_cortical_cell_types','var');
        limit_to_cortical_cell_types = false;
    end
    
    seperation_indices = length(giger_cell_types.samples_id);
    cell_type_order = {'Neurons'};
    
    
    [reorder_benchmark_cell_type, reorder_data_genes] = reorderUsingId(giger_cell_types.all_symbols, gene_info.gene_symbols);
    reorder_benchmark_cell_type_expresion = giger_cell_types.expression(reorder_benchmark_cell_type,:);
    reorder_cellmix_expresion = cell_mix.celltype_profile(reorder_data_genes,:);

    corr_matrix = corr(reorder_benchmark_cell_type_expresion, reorder_cellmix_expresion, 'type','spearman');
    draw_imagesc_with_seperation(corr_matrix,seperation_indices,cell_type_order);

    ax = gca;
    ax.XTick = 1:length(cell_mix.cell_types);
    ax.XTickLabel = cell_mix.cell_types;
    figure;
%     single_sample_scatter(reorder_mouse_cell_type_expresion(:,28), reorder_cellmix_expresion(:,strcmp('Neurons',cell_mix.cell_types)));
    
    sample_cell_type_id =double(giger_cell_types.sample2type) * ( ( 1:size(giger_cell_types.sample2type,2))');
    giger_cell_types.is_neuron = giger_cell_types.is_neuron(sample_cell_type_id);
    giger_cell_types.is_astro = giger_cell_types.is_astro(sample_cell_type_id);
    giger_cell_types.is_oligo = giger_cell_types.is_oligo(sample_cell_type_id);
%     aucs = drawROC(cell_mix, giger_cell_types, corr_matrix);
    [median_within, median_outside] = compute_median_corr(corr_matrix, giger_cell_types);
    
end

function [median_within, median_outside] = compute_median_corr(corr_matrix, giger_cell_types)
    astro_median_same = median( corr_matrix(giger_cell_types.is_astro,1));
    astro_median_diff = median( corr_matrix(~giger_cell_types.is_astro,1));
    
    neuron_median_same = median( corr_matrix(giger_cell_types.is_neuron,2));
    neuron_median_diff = median( corr_matrix(~giger_cell_types.is_neuron,2));
    
    oligo_median_same = median( corr_matrix(giger_cell_types.is_oligo,3));
    oligo_median_diff = median( corr_matrix(~giger_cell_types.is_oligo,3));
    
    median_within = [astro_median_same, neuron_median_same, oligo_median_same];
    median_outside = [astro_median_diff, neuron_median_diff, oligo_median_diff];
end
   
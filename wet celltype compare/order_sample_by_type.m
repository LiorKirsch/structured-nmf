function [wet_cell_types,seperation_indices,cell_type_order] = order_sample_by_type(wet_cell_types,limit_to_cortical_cell_types)
    [~,cell_type_ind] = max(double([wet_cell_types.is_neuron, wet_cell_types.is_astro, wet_cell_types.is_oligo]),[],2);
    [~,reorder_ind] = sort(cell_type_ind);
    
    [~,samples_cell_type_ind] = max(double([double(wet_cell_types.sample2type) * double(wet_cell_types.is_neuron), double(wet_cell_types.sample2type) * double(wet_cell_types.is_astro), double(wet_cell_types.sample2type) * double(wet_cell_types.is_oligo)]),[],2);
    [~,samples_reorder] = sort(samples_cell_type_ind);
%     samples_reorder = double(mouse_cell_types.sample2type) * reorder_ind;
%     [~, samples_reorder] = sort(samples_reorder);
    
    wet_cell_types = limit_data_by_cell_type_filter(wet_cell_types, reorder_ind);
    wet_cell_types = limit_data_by_sample_filter(wet_cell_types, samples_reorder);    
    
    if (limit_to_cortical_cell_types)
       wet_cell_types = limit_to_samples_from_cortex(wet_cell_types);
    end
        
    ind_neuron = max(find(double(wet_cell_types.sample2type) * double(wet_cell_types.is_neuron)));
    ind_astro = max(find(double(wet_cell_types.sample2type) * double(wet_cell_types.is_astro)));
    ind_oligo = max(find(double(wet_cell_types.sample2type) * double(wet_cell_types.is_oligo)));
    cell_type_order = {'Neurons','Astrocytes','Oligodendrocytes'};
    seperation_indices = [ind_neuron,ind_astro, ind_oligo];
    fprintf('neurons -> %d, astro -> %d, oligo -> %d\n',ind_neuron, ind_astro, ind_oligo);
end


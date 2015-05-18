function mouse_cell_types = limit_data_by_cell_type_filter(mouse_cell_types, cell_type_filter)
    mouse_cell_types.cell_type_description = mouse_cell_types.cell_type_description(cell_type_filter);
    mouse_cell_types.anatomical_region = mouse_cell_types.anatomical_region(cell_type_filter);
    mouse_cell_types.age_postnatal_day = mouse_cell_types.age_postnatal_day(cell_type_filter);
    mouse_cell_types.method_for_single_cell = mouse_cell_types.method_for_single_cell(cell_type_filter);
    mouse_cell_types.RNA_isolation_method = mouse_cell_types.RNA_isolation_method(cell_type_filter);
    mouse_cell_types.RNA_amplification_and_labeling_method = mouse_cell_types.RNA_amplification_and_labeling_method(cell_type_filter);
    mouse_cell_types.RNA_input_amount_to_microarray = mouse_cell_types.RNA_input_amount_to_microarray(cell_type_filter);
    mouse_cell_types.microarray_platform = mouse_cell_types.microarray_platform(cell_type_filter);
    mouse_cell_types.reference = mouse_cell_types.reference(cell_type_filter);
    mouse_cell_types.cell_type_id = mouse_cell_types.cell_type_id(cell_type_filter);
    
    mouse_cell_types.is_neuron = mouse_cell_types.is_neuron(cell_type_filter);
    mouse_cell_types.is_astro = mouse_cell_types.is_astro(cell_type_filter);
    mouse_cell_types.is_oligo = mouse_cell_types.is_oligo(cell_type_filter);
    mouse_cell_types.is_cortex_or_hippocampus = mouse_cell_types.is_cortex_or_hippocampus(cell_type_filter);
    
    mouse_cell_types.sample2type = mouse_cell_types.sample2type(:, cell_type_filter);
end

function mouse_cell_types = limit_data_by_sample_filter(mouse_cell_types, sample_filter)    
    mouse_cell_types.expression = mouse_cell_types.expression(:,sample_filter);
    mouse_cell_types.sample2type = mouse_cell_types.sample2type(sample_filter, :);
    mouse_cell_types.samples_id = mouse_cell_types.samples_id(sample_filter);
end

function [neuro_inds, oligo_inds, astro_inds] = get_celltype_inds(...
    mouse_cell_types, filters_string)

    num_types = length(mouse_cell_types.is_neuron);
    mask = true(num_types,1);
    
    filters = strsplit(filters_string,'_');
    for i = 1:length(filters)
        switch filters{i}
            case 'all'
                % do nothing
            case 'cortex' 
                cortex_inds = mouse_cell_types.is_cortex_or_hippocampus;
                mask = mask & cortex_inds;
            case 'doyle'
                doyle_inds = strmatch('Doyle',mouse_cell_types.reference);
                doyle_mask = false(num_types,1);
                doyle_mask(doyle_inds) = true;
                mask = mask & doyle_mask;
            otherwise
                error('filter not supported - %s', filters{i} );
        end
    end
    
    neuro_inds =  double(mouse_cell_types.is_neuron & mask);
    oligo_inds =  double(mouse_cell_types.is_oligo & mask);
    astro_inds =  double(mouse_cell_types.is_astro & mask);
end
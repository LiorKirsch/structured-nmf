function parms = get_structure_matrix(dataset, structure_type,selected_regions, parms)


switch dataset
    case {'kang_all_regions','Human6_selected_regions', 'brainspan_rnaseq'}
       
    case 'Zapala_isocortex_medulla_striatum_cerebellum'
        [parms.structure_matrix,parms.relation_regions] = zapala_tree_structure(selected_regions);
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);

    otherwise
        error('unkown dataset %s - cannot strucutre matrix',dataset)
end


end
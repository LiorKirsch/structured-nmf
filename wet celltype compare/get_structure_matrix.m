function parms = get_structure_matrix(dataset, structure_type,selected_regions, parms)

num_regions = length(selected_regions);

switch dataset
    case 'Human6_selected_regions'
        error('no tree for human 6 yet');
    case {'kang_all_regions','brainspan_rnaseq_DFC_OFC'}
        [parms.structure_matrix,parms.relation_regions] = kang_tree_structure(selected_regions);
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);
            
    case 'Zapala_isocortex_medulla_striatum_cerebellum'
        [parms.structure_matrix,parms.relation_regions] = zapala_tree_structure(selected_regions);
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);

    otherwise
        error('unkown dataset %s - cannot strucutre matrix',dataset)
end


end
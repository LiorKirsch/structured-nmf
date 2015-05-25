function parms = get_structure_matrix(dataset, structure_type,selected_regions, parms)

num_regions = length(selected_regions);

dataset = strsplit(dataset,'_');
dataset = dataset{1};
switch dataset
    case 'Human6'
        addpath('/home/lab/lior/Projects/buildStructureOntology/');
        selected_ontology = load('/home/lab/lior/Projects/buildStructureOntology/humanOntologyObject.mat', 'humanOntology');
        selected_ontology = selected_ontology.humanOntology;
        tree_structure_matrix = selected_ontology.dependencyMatrix;
        node_child_recursive =  inv(eye(size(tree_structure_matrix)) - tree_structure_matrix); % including self
        limit_to = ismember(selected_ontology.structureLabels(:,4), selected_regions);

        keep_strct = any(node_child_recursive(:,limit_to),2);
        region_with_parents = selected_ontology.structureLabels(keep_strct,4);
        tree_structure_matrix = tree_structure_matrix(keep_strct,keep_strct);
        
        parms.structure_matrix = tree_structure_matrix;
        parms.relation_regions = region_with_parents;
        
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);
        
    case {'kang','brainspan'}
        [parms.structure_matrix,parms.relation_regions] = kang_tree_structure(selected_regions);
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);
            
    case 'Zapala'
        [parms.structure_matrix,parms.relation_regions] = zapala_tree_structure(selected_regions);
        [parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,selected_regions,structure_type);

    otherwise
        error('unkown dataset %s - cannot strucutre matrix',dataset)
end


end
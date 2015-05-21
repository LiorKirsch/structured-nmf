function [tree_structure_matrix,all_strct] = zapala_tree_structure(limit_to)
%     Nervous_system -> Telencephalon
%     Nervous_system -> Diencephalon
%     Nervous_system -> Mesencephalon
%     Nervous_system -> Metencephalon
%     Nervous_system -> Myelencephalon
%     
%     Telencephalon -> Olfactory Bulbs
%     Telencephalon -> Isocortex

%     Isocortex -> Motor Cortex
%     Isocortex -> Entorhinal Cortex
%     Isocortex -> Perirhinal Cortex
%     Isocortex -> Hippocampus

%     Hippocampus -> Hippocampus CA1
%     Hippocampus -> Hippocampus CA3
%     Hippocampus -> Dentate Gyrus

%     Telencephalon -> Cerebral Nuclei
%     Cerebral Nuclei -> Amygdala
%     Cerebral Nuclei -> Striatum

%     Diencephalon -> Thalamus
%     Diencephalon -> Hypothalamus
%     Diencephalon -> Pituitary

%     Mesencephalon -> Inferior Colliculi
%     Mesencephalon -> Superior Colliculi

%     Metencephalon -> cerebellum
%     Metencephalon -> Pons

%     Myelencephalon -> Medulla
    
    
    all_strct = {'Nervous_system','Telencephalon','Diencephalon','Mesencephalon',...
        'Metencephalon','Myelencephalon','Cerebral_cortex','Olfactory Bulbs',...
        'Motor Cortex','Entorhinal Cortex','Perirhinal Cortex','Hippocampus','Hippocampus CA1'...
        'Hippocampus CA3','Dentate Gyrus','Cerebral Nuclei','Amygdala','Striatum',...
        'Thalamus','Hypothalamus','Pituitary','Inferior Colliculi','Superior Colliculi',...
        'Cerebellum','Pons','Medulla'}';

    
        % add BNST 'Bed nucleus of the stria terminalis'
    
    tree_structure_matrix = false(length(all_strct),length(all_strct));
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'Telencephalon')) = true;
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'Diencephalon')) = true;
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'Mesencephalon')) = true;
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'Metencephalon')) = true;
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'Myelencephalon')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'Olfactory Bulbs')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'Cerebral_cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'Motor Cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'Entorhinal Cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'Perirhinal Cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'Hippocampus')) = true;
    tree_structure_matrix(strcmp(all_strct,'Hippocampus'),...
                          strcmp(all_strct,'Hippocampus CA1')) = true;
    tree_structure_matrix(strcmp(all_strct,'Hippocampus'),...
                          strcmp(all_strct,'Hippocampus CA3')) = true;
    tree_structure_matrix(strcmp(all_strct,'Hippocampus'),...
                          strcmp(all_strct,'Dentate Gyrus')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'Cerebral Nuclei')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral Nuclei'),...
                          strcmp(all_strct,'Amygdala')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral Nuclei'),...
                          strcmp(all_strct,'Striatum')) = true;
    tree_structure_matrix(strcmp(all_strct,'Diencephalon'),...
                          strcmp(all_strct,'Thalamus')) = true;
    tree_structure_matrix(strcmp(all_strct,'Diencephalon'),...
                          strcmp(all_strct,'Hypothalamus')) = true;
    tree_structure_matrix(strcmp(all_strct,'Diencephalon'),...
                          strcmp(all_strct,'Pituitary')) = true;
    tree_structure_matrix(strcmp(all_strct,'Mesencephalon'),...
                          strcmp(all_strct,'Inferior Colliculi')) = true;
    tree_structure_matrix(strcmp(all_strct,'Mesencephalon'),...
                          strcmp(all_strct,'Superior Colliculi')) = true;
    tree_structure_matrix(strcmp(all_strct,'Metencephalon'),...
                          strcmp(all_strct,'Cerebellum')) = true;
    tree_structure_matrix(strcmp(all_strct,'Metencephalon'),...
                          strcmp(all_strct,'Pons')) = true;
    tree_structure_matrix(strcmp(all_strct,'Myelencephalon'),...
                          strcmp(all_strct,'Medulla')) = true;

    if exist('limit_to','var')
        node_child_recursive =  inv(eye(size(tree_structure_matrix)) - tree_structure_matrix); % including self
        limit_to = ismember(all_strct, limit_to);

        keep_strct = any(node_child_recursive(:,limit_to),2);
        all_strct = all_strct(keep_strct);
        tree_structure_matrix = tree_structure_matrix(keep_strct,keep_strct);
    end
    
end

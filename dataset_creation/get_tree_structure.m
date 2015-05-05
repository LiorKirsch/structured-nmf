function [tree_structure_matrix,all_strct] = get_tree_structure(shortcut_only_child,limit_to)
%     Nervous_system -> Telencephalon
%     Nervous_system -> Diencephalon
%     Nervous_system -> Mesencephalon
%     Nervous_system -> Metencephalon
%     Nervous_system -> Myelencephalon
%     Nervous_system -> spinal_cord
%     
%     Telencephalon -> Cerebral_cortex
%     Cerebral_cortex -> cortex_L5A
%     Cerebral_cortex -> cortex_L5B
%     Cerebral_cortex -> cortex_L6
%     Telencephalon -> striatum
%     Metencephalon -> cerebellum
%     Myelencephalon -> brainstem
    
    
    all_strct = {'Nervous_system','Telencephalon','Diencephalon','Mesencephalon',...
        'Metencephalon','Myelencephalon','spinal_cord','Cerebral_cortex','cortex_L5A',...
        'cortex_L5B','cortex_L6','striatum','cerebellum','brainstem'}';

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
    tree_structure_matrix(strcmp(all_strct,'Nervous_system'),...
                          strcmp(all_strct,'spinal_cord')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'Cerebral_cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'cortex_L5A')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'cortex_L5B')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral_cortex'),...
                          strcmp(all_strct,'cortex_L6')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'striatum')) = true;
    tree_structure_matrix(strcmp(all_strct,'Metencephalon'),...
                          strcmp(all_strct,'cerebellum')) = true;
    tree_structure_matrix(strcmp(all_strct,'Myelencephalon'),...
                          strcmp(all_strct,'brainstem')) = true;

    if exist('limit_to','var')
        node_child_recursive =  inv(eye(size(tree_structure_matrix)) - tree_structure_matrix); % including self
        limit_to = ismember(all_strct, limit_to);

        keep_strct = any(node_child_recursive(:,limit_to),2);
        all_strct = all_strct(keep_strct);
        tree_structure_matrix = tree_structure_matrix(keep_strct,keep_strct);
    end
    
    if ~exist('shortcut_only_child','var')
        shortcut_only_child = false;
    end
    
    if shortcut_only_child
        only_child = find(sum(tree_structure_matrix,2) ==1);
        for i =1:length(only_child)
            curr_ind = only_child(i);
            child_ind = find(tree_structure_matrix( curr_ind,:));
            parent_ind = find(tree_structure_matrix( :,curr_ind));
            tree_structure_matrix(parent_ind, child_ind) = true;
        end
        tree_structure_matrix(only_child,:) = [];
        tree_structure_matrix(:,only_child) = [];
        all_strct(only_child) = [];
    end
end

function [tree_structure_matrix,all_strct] = kang_tree_structure(limit_to)
%     Nervous_system -> Telencephalon
%     Nervous_system -> Diencephalon
%     Nervous_system -> Mesencephalon
%     Nervous_system -> Metencephalon
%     Nervous_system -> Myelencephalon
%     
%     Telencephalon -> cerebellar cortex
%     Telencephalon -> Cerebral Nuclei

%     cerebellar cortex -> primary auditory (A1) cortex
%     cerebellar cortex -> posterior inferior parietal cortex
%     cerebellar cortex -> primary motor (M1) cortex
%     cerebellar cortex -> hippocampus
%     cerebellar cortex -> primary somatosensory (S1) cortex
%     cerebellar cortex -> primary visual (V1) cortex
%     cerebellar cortex -> prefrontal cortex
%     cerebellar cortex -> temporal cortex'

%     temporal cortex -> superior temporal cortex
%     temporal cortex -> inferior temporal cortex

%     prefrontal cortex -> medial prefrontal cortex 
%     prefrontal cortex -> orbital prefrontal cortex
%     prefrontal cortex -> ventrolateral prefrontal cortex
%     prefrontal cortex -> dorsolateral prefrontal cortex

%     Cerebral Nuclei -> amygdala
%     Cerebral Nuclei -> striatum

%     Diencephalon -> 'mediodorsal nucleus of the thalamus'

    
    
    all_strct = {'Nervous_system','Telencephalon','Diencephalon','Mesencephalon',...
        'Metencephalon','cerebellar cortex','primary auditory (A1) cortex',...
    'hippocampus','posterior inferior parietal cortex',...
    'primary motor (M1) cortex',...
    'primary somatosensory (S1) cortex', 'primary visual (V1) cortex',...
    'temporal cortex', 'superior temporal cortex','inferior temporal cortex',...
    'prefrontal cortex', 'medial prefrontal cortex','orbital prefrontal cortex',...
    'ventrolateral prefrontal cortex','dorsolateral prefrontal cortex',...
    'Cerebral Nuclei','striatum','amygdala',...
    'mediodorsal nucleus of the thalamus'};

    
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
                          strcmp(all_strct,'cerebellar cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'Telencephalon'),...
                          strcmp(all_strct,'Cerebral Nuclei')) = true;
                      
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'primary auditory (A1) cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'posterior inferior parietal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'primary motor (M1) cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'hippocampus')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'primary somatosensory (S1) cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'primary visual (V1) cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'temporal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'cerebellar cortex'),...
                          strcmp(all_strct,'prefrontal cortex')) = true;
    
    tree_structure_matrix(strcmp(all_strct,'temporal cortex'),...
                          strcmp(all_strct,'superior temporal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'temporal cortex'),...
                          strcmp(all_strct,'inferior temporal cortex')) = true;
                      
    tree_structure_matrix(strcmp(all_strct,'prefrontal cortex'),...
                          strcmp(all_strct,'medial prefrontal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'prefrontal cortex'),...
                          strcmp(all_strct,'orbital prefrontal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'prefrontal cortex'),...
                          strcmp(all_strct,'ventrolateral prefrontal cortex')) = true;
    tree_structure_matrix(strcmp(all_strct,'prefrontal cortex'),...
                          strcmp(all_strct,'dorsolateral prefrontal cortex')) = true;
   
    tree_structure_matrix(strcmp(all_strct,'Cerebral Nuclei'),...
                          strcmp(all_strct,'amygdala')) = true;
    tree_structure_matrix(strcmp(all_strct,'Cerebral Nuclei'),...
                          strcmp(all_strct,'striatum')) = true;
    tree_structure_matrix(strcmp(all_strct,'Diencephalon'),...
                          strcmp(all_strct,'mediodorsal nucleus of the thalamus')) = true;
   

    if exist('limit_to','var')
        node_child_recursive =  inv(eye(size(tree_structure_matrix)) - tree_structure_matrix); % including self
        limit_to = ismember(all_strct, limit_to);

        keep_strct = any(node_child_recursive(:,limit_to),2);
        all_strct = all_strct(keep_strct);
        tree_structure_matrix = tree_structure_matrix(keep_strct,keep_strct);
    end
    
end


function [gene_inds_target, gene_inds_mouse] = ...
        get_intersecting_genes(gene_symbols_target, gene_symbols_mouse, parms)
% 
    switch parms.species
        case 'mouse'
            [~, gene_inds_target, gene_inds_mouse] = ...
                intersect(gene_symbols_target, gene_symbols_mouse);

        case {'monkey', 'human'}
            addpath('/cortex/code/matlab/homologous_gene_mapping/');
            switch parms.species
                case 'monkey'
                    [gene_to_group_mouse, gene_to_group_target, ...
                     homologous_group_id] = ...
                        gene_to_homolog_group('mouse_laboratory', ...
                                              'rhesus_macaque', ...
                                              gene_symbols_mouse, ...
                                              'symbol', ...
                                              gene_symbols_target, 'symbol');
                case 'human'
                    [gene_to_group_mouse, gene_to_group_target, ...
                     homologous_group_id] = ...
                        gene_to_homolog_group('mouse_laboratory','human', ...
                                              gene_symbols_mouse, ...
                                              'symbol', ...
                                              gene_symbols_target,'symbol');
            end
            
%             groups_with_1_to_1 = sum(gene_to_group_mouse,1) == 1  & sum(gene_to_group_target,1) == 1;
%             gene_to_group_mouse = gene_to_group_mouse(:,groups_with_1_to_1);
%             gene_to_group_target = gene_to_group_target(:,groups_with_1_to_1);
%             gene_to_group_mouse = (1:size(gene_to_group_mouse,1)) * gene_to_group_mouse ;
%             gene_to_group_target = (1:size(gene_to_group_target,1)) * gene_to_group_target ;

%             gene_inds_mouse = gene_to_group_mouse;
%             gene_inds_target = gene_to_group_target;
            
            gene_inds_mouse = nan;
            gene_inds_target = any(gene_to_group_target,2);


    end
    
end
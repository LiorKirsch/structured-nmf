function [gene_symbol_output,scores_output] = get_homolog_symbols_and_scores(gene_symbols_mouse,gene_scores, parms)
% 

            addpath('/cortex/code/matlab/homologous_gene_mapping/');
            
            [gene_to_group_mouse, gene_to_group_target, ...
             homologous_group_id,gene_id_mouse,gene_id_target] = ...
                gene_to_homolog_group('mouse_laboratory','human', ...
                                      gene_symbols_mouse, ...
                                      'symbol');
            
            num_genes_in_group = sum(gene_to_group_mouse,1);
            avg_homologus_group_score = gene_scores * gene_to_group_mouse;
            avg_homologus_group_score = avg_homologus_group_score ./ num_genes_in_group;
            
            
            scores_output = gene_to_group_target * avg_homologus_group_score';
            gene_symbol_output = gene_id_target;

    
end
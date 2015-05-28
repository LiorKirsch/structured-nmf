function top_gene_symbols = okaty_select_genes2(num_top_genes, ...
                                               selection_score, parms)
%
% Select genes basedon Barres2014. 
%
   
        gene_okaty_filter = take_from_struct(parms, 'gene_okaty_filter', 'all');
        pattern = '([a-z_]*)(\d*)';
        [tokens, match] = regexp(parms.gene_subset, pattern, 'tokens', 'match');
        selection_method = tokens{1}{1};    
        curr_parms.gene_subset = selection_method;
        curr_parms.gene_okaty_filter = gene_okaty_filter;
        gene_subset_file = set_filenames('gene_subset', curr_parms);
    
        vars = {'sorted_symbols_human','sorted_scores_human',...
                'sorted_symbols_mouse','sorted_scores_mouse'};
        if exist(gene_subset_file,'file')
            fprintf('loading info gain from disk - %s\n', gene_subset_file);
            load(gene_subset_file, vars{1:end});
        else
            mouse_cell_types = load('mouse_cell_type_profiles.mat');
            [neuro_inds, oligo_inds, astro_inds] = get_celltype_inds(...
                mouse_cell_types, gene_okaty_filter);
            
            sample_neuro = mouse_cell_types.sample2type * neuro_inds;
            sample_oligo = mouse_cell_types.sample2type * oligo_inds;
            sample_astro = mouse_cell_types.sample2type * astro_inds;            
            expression = mouse_cell_types.expression';

            y = sample_neuro * 1 + sample_astro* 2 + sample_oligo *3;
            keep_inds = y>0;
            y = y(keep_inds);
            expression = expression(keep_inds,:);
            
            
            switch selection_score
              case 'anova',
                for i=1:size(expression,2)
                    printPercentCounter(i, size(expression,2));
                    p(i) = anova1(expression(:,i), y, 'off');
                end
                sort_direction = 'ascend';
                [sorted_p, sort_inds] = sort(p);
                [sorted_scores_mouse, sort_all_names_inds] = sort(p(mouse_cell_types.refer_to_index),sort_direction);
              case  'infogain',
                sort_direction = 'descend';
                y = arrayfun(@(x) sprintf('%d',x),y,'Uniformoutput',false);
                [~, gainAttrs, ~, sort_inds] = infoGain(expression, y);
                [sorted_scores_mouse, sort_all_names_inds] = sort( gainAttrs(...
                mouse_cell_types.refer_to_index),sort_direction);
              case  'gainratio',
                sort_direction = 'descend';
                y = arrayfun(@(x) sprintf('%d',x),y,'Uniformoutput',false);
                
                [igsp, iga, sig, siga] = infoGain( expression,y );
        [ ~,gainAttrs,~,sort_inds ] = gainRatio(expression,y,iga);
                [sorted_scores_mouse, sort_all_names_inds] = sort( gainAttrs(...
                mouse_cell_types.refer_to_index),sort_direction);
              otherwise, 
                error('invalid selection_score = [%s]\n', selection_score)
            end
%             sorted_symbols = mouse_cell_types.gene_symbol(sort_inds);
            sorted_symbols_mouse = mouse_cell_types.all_symbols(sort_all_names_inds);

            [symbols_human,scores_human] = get_homolog_symbols_and_scores(sorted_symbols_mouse,sorted_scores_mouse, parms);
            
            [sorted_scores_human, sort_scores_human_inds] = sort( scores_human,sort_direction);
            sorted_symbols_human = symbols_human(sort_scores_human_inds);
            save(gene_subset_file, vars{1:end});
        end
        
    
        switch parms.species
            case 'mouse'
                top_gene_symbols = sorted_symbols_mouse(1:num_top_genes);   
            case 'human'
                top_gene_symbols = sorted_symbols_human(1:num_top_genes);   
            otherwise
                error('species not supported - %s', parms.species);
        end
end


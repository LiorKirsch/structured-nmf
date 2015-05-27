function top_gene_symbols = okaty_select_genes(num_top_genes, ...
                                               selection_score, parms)
%
% Select genes basedon Barres2014. 
%
    persistent local_sorted_symbols
    if isempty(local_sorted_symbols)
    
        gene_okaty_filter = take_from_struct(parms, 'gene_okaty_filter', 'all');
        pattern = '([a-z_]*)(\d*)';
        [tokens, match] = regexp(parms.gene_subset, pattern, 'tokens', 'match');
        selection_method = tokens{1}{1};    
        curr_parms.gene_subset = selection_method;
        curr_parms.gene_okaty_filter = gene_okaty_filter;
        gene_subset_file = set_filenames('gene_subset', curr_parms)
    
        vars = {'sort_inds', 'mouse_cell_types', 'sorted_symbols'};
        if exist(gene_subset_file,'file')
            fprintf('load gene scores from [%s]\n', gene_subset_file);
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
                    if mod(i,100)==0,  fprintf('.'); end
                    p(i) = anova1(expression(:,i), y, 'off');
                end
                [sorted_p, sort_inds] = sort(p);
              case  'infogain', 
                y = arrayfun(@(x) sprintf('%d',x),y,'Uniformoutput',false);
                [~, ~, ~, sort_inds] = infoGain(expression, y);
              otherwise, 
                error('invalid selection_score = [%s]\n', selection_score)
            end
            sorted_symbols = mouse_cell_types.gene_symbol(sort_inds);
            save(gene_subset_file, vars{1:end});
        end
        local_sorted_symbols = sorted_symbols;
    else
        sorted_symbols = local_sorted_symbols;
    end
   top_gene_symbols = sorted_symbols(1:num_top_genes);   
end


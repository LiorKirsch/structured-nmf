function top_gene_symbols = okaty_infogain_genes(num_top_genes, parms)
%
% Select genes basedon Barres2014. 
%

   persistent local_sort_inds
   persistent local_mouse_cell_types
   
   if isempty(local_sort_inds)
        mouse_cell_types = load('mouse_cell_type_profiles.mat');
        sample_neuro = mouse_cell_types.sample2type * double(mouse_cell_types.is_neuron);
        sample_oligo = mouse_cell_types.sample2type * double(mouse_cell_types.is_oligo);
        sample_astro = mouse_cell_types.sample2type * double(mouse_cell_types.is_astro);

        expression = mouse_cell_types.expression';

        y = sample_neuro * 1 + sample_astro* 2 + sample_oligo *3;
        remove_inds = y==0;
        y = y(~remove_inds);
        expression = expression(~remove_inds,:);
        y = arrayfun(@(x) sprintf('%d',x),y,'Uniformoutput',false);

        [ ~,~,~,sortedAttr ] = infoGain( expression,y );
      
      local_sort_inds = sortedAttr;
      local_mouse_cell_types = mouse_cell_types;
   end

   top_gene_symbols = local_mouse_cell_types.gene_symbol(local_sort_inds(1:num_top_genes));
   
end
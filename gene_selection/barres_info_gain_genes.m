function top_gene_symbols = barres_info_gain_genes(num_top_genes, parms)
%
% Select genes basedon Barres2014. 
%

   persistent local_sort_inds
   persistent local_barres
   
   if isempty(local_sort_inds)

       % load Barres data 
       filename  = fullfile('/','cortex','data','RNA-Seq','mouse', ...
                           'Barres-2014','barres_rnaseq.mat');
       barres = load(filename);
       
       [ ~,~,~,sortedAttr ] = infoGain( barres.data',barres.cell_types' );
       local_sort_inds = sortedAttr;
       local_barres = barres;
   end

   top_gene_symbols = local_barres.gene_symbols(local_sort_inds(1:num_top_genes));
   
end
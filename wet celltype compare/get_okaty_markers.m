function [neuronal_markers,astro_markers,oligo_markers] = get_okaty_markers(k_markers, treshold)

 load('mouse_cell_type_profiles.mat');

expression = 2 .^ expression;



neuro_avg_exp = expression * sample2type * double(is_neuron) / sum(is_neuron);
astro_avg_exp = expression * sample2type * double(is_astro) / sum(is_astro);
oligo_avg_exp = expression * sample2type * double(is_oligo) / sum(is_oligo);
tmp = [neuro_avg_exp, astro_avg_exp, oligo_avg_exp];



small_inds = tmp(:,2) < treshold &    tmp(:,3) < treshold ;
small_genes = tmp(small_inds,:);
tmp_symbols = gene_symbol(small_inds);

[a,a_sort_inds] = sort( small_genes(:,1),'descend');
neuronal_markers = tmp_symbols(a_sort_inds(1:k_markers));


small_inds = tmp(:,1) < treshold &    tmp(:,3) < treshold ;
small_genes = tmp(small_inds,:);
tmp_symbols = gene_symbol(small_inds);

[a,a_sort_inds] = sort( small_genes(:,2),'descend');
astro_markers = tmp_symbols(a_sort_inds(1:k_markers));


small_inds = tmp(:,1) < treshold &    tmp(:,2) < treshold ;
small_genes = tmp(small_inds,:);
tmp_symbols = gene_symbol(small_inds);

[a,a_sort_inds] = sort( small_genes(:,3),'descend');
oligo_markers = tmp_symbols(a_sort_inds(1:k_markers));


% get the correct names - for genes which have duplicate names in
% gene_symbol
neuronal_markers = all_symbols(ismember(refer_to_index, find( ismember(gene_symbol, neuronal_markers) )));
astro_markers = all_symbols(ismember(refer_to_index, find( ismember(gene_symbol, astro_markers) )));
oligo_markers = all_symbols(ismember(refer_to_index, find( ismember(gene_symbol, oligo_markers) )));

end
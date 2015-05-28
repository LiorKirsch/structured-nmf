function H_markers = get_markers_from_GT(GT_profiles, k_markers, treshold,parms)

num_genes = size(GT_profiles,1);

small_inds = find( GT_profiles(:,2) < treshold &    GT_profiles(:,3) < treshold );
small_genes = GT_profiles(small_inds,:);

[~,a_sort_inds] = sort( small_genes(:,1),'descend');
neuronal_markers_inds = small_inds(a_sort_inds(1:k_markers));


small_inds = find(GT_profiles(:,1) < treshold &    GT_profiles(:,3) < treshold) ;
small_genes = GT_profiles(small_inds,:);

[a,a_sort_inds] = sort( small_genes(:,2),'descend');
astro_markers_inds = small_inds(a_sort_inds(1:k_markers));


small_inds = find(GT_profiles(:,1) < treshold &    GT_profiles(:,2) < treshold) ;
small_genes = GT_profiles(small_inds,:);


[a,a_sort_inds] = sort( small_genes(:,3),'descend');
oligo_markers_inds = small_inds(a_sort_inds(1:k_markers));


H_markers = false(parms.num_types,num_genes);
H_markers(1, neuronal_markers_inds ) = true;
H_markers(2, astro_markers_inds ) = true;
H_markers(3, oligo_markers_inds ) = true;

non_unique = sum(H_markers,1) > 1;
H_markers(:,non_unique) = false;

end
function [gene_info, expression, parms] = gene_subset_selection(gene_info, ...
                                                      expression, parms)
%
    gene_subset = take_from_struct(parms, 'gene_subset', 'barres100')
    num_genes = length(gene_info.gene_symbols);
    switch gene_subset

        case 'all', 
            gene_mask = true(num_genes,1);
        
        case 'allen_subset', % 'genes_with_orthologs'
            allen_mouse_genes = load('allen_mouse_genes');
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbolst, allen_mouse_genes, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;

        case 'orth', % 'genes_with_orthologs'
            mouse_cell_types = load('mouse_cell_type_profiles.mat');
            species = take_from_struct(parms, 'species');    
            [~,gene_inds] = compare_to_true_profile(mouse_cell_types, ...
                                                    gene_info, species, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;

        case 'barres_discrim100' % top 100 discriminative barres
            top_gene_symbols = load_top_barres_genes(100, parms);
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;

        case 'barres_discrim1000' % top 1000 discriminative barres
            top_gene_symbols = load_top_barres_genes(1000, parms);
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
        
        case 'okaty_discrim100' % top 1000 discriminative barres
            top_gene_symbols = load_top_okaty_genes(100, parms);
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
            
        case 'okaty_discrim1000' % top 1000 discriminative barres
            top_gene_symbols = load_top_okaty_genes(1000, parms);
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
        
        case 'okaty_infogain_1000' % top 1000 infogain okaty
            top_gene_symbols = okaty_infogain_genes(1000, parms);
%             [~, ~, gene_inds] = intersect(upper(top_gene_symbols), gene_info.gene_symbols);
             gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
         
         case 'okaty_gainratio_1000' % top 1000 okaty
            top_gene_symbols = okaty_gainratio_genes(1000, parms);
%             [~, ~, gene_inds] = intersect(upper(top_gene_symbols), gene_info.gene_symbols);
            gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
           
         case 'barres_infogain_1000' % top 1000 infogain barres
            top_gene_symbols = barres_info_gain_genes(1000, parms);
%             [~, ~, gene_inds] = intersect(upper(top_gene_symbols), gene_info.gene_symbols);
             gene_inds = get_intersecting_genes(...
                gene_info.gene_symbols, top_gene_symbols, parms);
            gene_mask = false(num_genes, 1);
            gene_mask(gene_inds) = true;
        
        otherwise
            error('unkown gene_subset = [%s]', gene_subset);
    end

    gene_info.gene_symbols = gene_info.gene_symbols(gene_mask);
    gene_info.entrez_ids = gene_info.entrez_ids(gene_mask);    
    expression = expression(:,gene_mask);
    parms.gene_hash = sum( gene_mask .* ((1:num_genes)') );
end

% ========================================================
function top_gene_symbols = load_top_barres_genes(num_top_genes, parms)
%
% Select genes basedon Barres2014. 
%

   persistent local_top_gene_symbols
   persistent local_num_top_genes
   if isempty(local_top_gene_symbols) || ...
           isempty(intersect(local_num_top_genes,num_top_genes))

       % load Barres data 
       filename  = fullfile('/','cortex','data','RNA-Seq','mouse', ...
                           'Barres-2014','barres_rnaseq.mat');
       barres = load(filename);
       inds = discriminative_feature_score(barres.data);
       
       local_top_gene_symbols = barres.gene_symbols(inds(1:num_top_genes));
       local_num_top_genes = num_top_genes;
   end

   top_gene_symbols = local_top_gene_symbols;

end
function top_inds = discriminative_feature_score(data)
       [num_genes, num_tissues] = size(data);
       p = data/sum(data(:));
       mi = zeros(num_genes, 1);

       py = sum(p,1);    
       hy = -sum(py.*log(py));
       fprintf('    ');
       for i_gene = 1:num_genes
           printPercentCounter(i_gene, num_genes);
           pp = nan(2, num_tissues);
           pp(1,:) = p(i_gene, :);
           pp(2,:) = 1-p(i_gene, :);
           mi(i_gene) = compute_mi(pp, hy);
       end
       [sorted, top_inds] = sort(mi, 'descend');
end
function mi = compute_mi(pxy, hy)
    px = sum(pxy,2);
    p1 = pxy(1,:) / sum(pxy(1,:));
    p2 = pxy(2,:) / sum(pxy(2,:));

    hy_given_x(1) = -sum(p1.*log(p1));
    hy_given_x(2) = -sum(p2.*log(p2));

    %py = sum(pxy,1);    
    %hy = -sum(py.*log(py))
    mi = hy - hy_given_x*px;
end


% ========================================================
function top_gene_symbols = load_top_okaty_genes(num_top_genes, parms)
%
% Select genes basedon Okaty2014. 
%

   persistent local_top_gene_symbols
   persistent local_num_top_genes
   if isempty(local_top_gene_symbols) || ...
           isempty(intersect(local_num_top_genes,num_top_genes))

       % load Okaty data 
       okaty_celltypes = load('mouse_cell_type_profiles.mat');
       mean_expression = okaty_celltypes.expression * okaty_celltypes.sample2type ;
       mean_expression = bsxfun(@rdivide, mean_expression, sum( okaty_celltypes.sample2type,1));
       
       inds = discriminative_feature_score(mean_expression);
       
       local_top_gene_symbols = okaty_celltypes.gene_symbol(inds(1:num_top_genes));
       local_num_top_genes = num_top_genes;
        
   end
   
   top_gene_symbols = local_top_gene_symbols;
end


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
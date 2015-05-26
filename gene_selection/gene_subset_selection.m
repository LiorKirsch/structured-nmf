function [gene_info, expression, parms] = gene_subset_selection(gene_info, ...
                                                      expression, parms)
%
    gene_subset = take_from_struct(parms, 'gene_subset', 'barres100');
%     gene_okaty_filter = take_from_struct(parms, 'gene_okaty_filter', 'all_types');
    num_genes = length(gene_info.gene_symbols);

    pattern = '([a-z_]*)(\d*)';
    [tokens, match] = regexp(gene_subset, pattern, 'tokens', 'match');
    selection_method = tokens{1}{1};
    if ~isempty(tokens{1}{2})
        num_to_select = str2num(tokens{1}{2});
    end
    
    
       switch selection_method
            case 'all', 
                gene_inds = (1:num_genes);

            case 'barres_discrim'
                top_gene_symbols = load_top_barres_genes(num_to_select, parms);

            case 'okaty_discrim'
                top_gene_symbols = load_top_okaty_genes(num_to_select, parms);

            case 'okaty_infogain'
                top_gene_symbols = okaty_infogain_genes(num_to_select, parms);

            case 'okaty_gainratio'
                top_gene_symbols = okaty_gainratio_genes(num_to_select, parms);

            case 'barres_infogain'
                top_gene_symbols = barres_info_gain_genes(1000, parms);

            case 'allen_subset', % ?
                allen_mouse_genes = load('allen_mouse_genes');
                gene_inds = get_intersecting_genes(...
                    gene_info.gene_symbols, allen_mouse_genes, parms);

            case 'orth', % 'genes_with_orthologs'
                mouse_cell_types = load('mouse_cell_type_profiles.mat');
                species = take_from_struct(parms, 'species');    
                [~,gene_inds] = compare_to_true_profile(mouse_cell_types, ...
                                                        gene_info, species, parms);
            otherwise
                error('unkown selection_method [%s], gene_subset = [%s]', ...
                      selection_method, gene_subset);
        end
        
       

    if ~exist('gene_inds', 'var')
        gene_inds = get_intersecting_genes(...
            gene_info.gene_symbols, top_gene_symbols, parms);
    end
        
    gene_info.gene_symbols = gene_info.gene_symbols(gene_inds);
    gene_info.entrez_ids = gene_info.entrez_ids(gene_inds);    
    expression = expression(:, gene_inds);
    parms.gene_hash = sum(gene_inds);
end

% ========================================================
function top_gene_symbols = load_top_barres_genes(num_top_genes, parms)
%
% Select genes basedon Barres2014. 
%
  
       % load Barres data 
       filename  = fullfile('/','cortex','data','RNA-Seq','mouse', ...
                           'Barres-2014','barres_rnaseq.mat');
       barres = load(filename);
       inds = discriminative_feature_score(barres.data);
       
       top_gene_symbols = barres.gene_symbols(inds(1:num_top_genes));
    
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

 
       % load Okaty data 
       okaty_celltypes = load('mouse_cell_type_profiles.mat');
       mean_expression = okaty_celltypes.expression * okaty_celltypes.sample2type ;
       mean_expression = bsxfun(@rdivide, mean_expression, sum( okaty_celltypes.sample2type,1));
       
       inds = discriminative_feature_score(mean_expression);
       
      top_gene_symbols = okaty_celltypes.gene_symbol(inds(1:num_top_genes));
    
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
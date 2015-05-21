function gene_entropy = find_genes_with_most_entropy(expression_matrix)

    [num_samples, num_genes] = size(expression_matrix);

    gene_entropy = nan(num_genes,1);
    for i = 1:num_genes
    gene_entropy(i) = wentropy(expression_matrix(:,i),'shannon');
    end

end

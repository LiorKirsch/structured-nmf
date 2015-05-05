function err = nmf_euclidean_dist(X,Y)
    err = sum(sum((X-Y).^2));
end
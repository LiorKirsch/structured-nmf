function H_marker_part = update_H_markers(X, W, H_markers,H_lambda,H_regularizer)
% H_markers is matrix which specifies for each type (K) which features are markers.
%
% TODO:   replace the "for loop" with matrix notations.
%
    K = size(H_markers,1);
    M = size(X,2);
    H_marker_part = zeros(K,M);
    WT_X = W'*X; % output size M
    sum_W = sum(W,1);  % output size K
    sum_W_square = sum(W.^2,1);  % output size K
    for i =1:K
        if H_lambda > 0
            H_marker_part(i,H_markers(i,:)) = WT_X(i,H_markers(i,:)) + H_lambda *H_regularizer(i,H_markers(i,:)) ;
            H_marker_part(i,H_markers(i,:)) = H_marker_part(i,H_markers(i,:)) / (sum_W_square(i) + H_lambda);
        else
            H_marker_part(i,H_markers(i,:)) =...
                 WT_X(i,H_markers(i,:))   / sum_W_square(i) ;
        end
    end

end
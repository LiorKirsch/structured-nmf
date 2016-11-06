function kkt_dist = get_stopping_cond(H,gradH, projection_type)


    switch projection_type
        case {'on_simplex' ,'on_simplex_with_noise'}
%             min_grad = min(gradH,[],1);
%             kkt_resdiue = H .* bsxfun(@minus, gradH, min_grad);
%             kkt_dist = norm(kkt_resdiue(:));

            alpha =  -1 * sum(H.*gradH,1);
            grad_plus_alpha = bsxfun(@plus, gradH, alpha);
            kkt_dist = norm(grad_plus_alpha( grad_plus_alpha < 0 ...
                | H >0));
            
%             alpha = bsxfun(@rdivide,sum(-gradH .*(H > 0)), sum(H>0) );
%             grad_plus_alpha = bsxfun(@plus, gradH, alpha);
%             kkt_dist = norm( grad_plus_alpha(H >0 | grad_plus_alpha<0) ) ;
                       
        case 'inside_simplex'
%             min_grad = min(gradH,[],1);
% %             min_grad = max(min_grad,0);
%             kkt_positive = H .* bsxfun(@minus, gradH, min_grad);
%             
%             dist_from_simplex  = sum(H,1) - 1;
%             kkt_simplex = dist_from_simplex .* min_grad;
%             
%             kkt_dist = sqrt(norm(kkt_positive(:))^2 + norm(kkt_simplex)^2);
            
            active_col = sum(H,1) < 1;
            
            kkt_dist_active = get_stopping_cond(H(:,active_col),gradH(:,active_col), 'positive');
            kkt_dist_not_act = get_stopping_cond(H(:,~active_col),gradH(:,~active_col), 'on_simplex');
            kkt_dist = sqrt(kkt_dist_active^2 + kkt_dist_not_act^2);
        case 'positive'
            kkt_dist = norm(gradH(gradH < 0 | H >0));
        otherwise
            error('unknown option for projection type - %s', projection_type);
    end
    
end
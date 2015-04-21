function [W, H, best_score, proportions_score] = match_profiles_to_gt(W, H, GT, GT_proportions, corr_type)
%
% Find a permutation of the rows of H that matches best the
% ground-truth data GT
%
% Use the Hungarian algorithem to solve the perfect matching problem
% Since we are interested in the maximum weight matching problem
%  (and not the minimum) we simply use -1 times the weights (the
%  correlation matrix).
%
   num_types = size(H, 1);  
   num_types_GT = size(GT, 1); 

   if ~exist('corr_type','var')
       corr_type = 'Pearson';
   end
   r = corr(H', GT','type',corr_type);
   
   if (num_types  < num_types_GT )
       all_perms = allcomb(repmat({1:num_types},num_types_GT,1));
       num_perms = size(all_perms,1);   
       scores = zeros(num_perms,1);
       for i=1:num_perms
           perm = all_perms(i,:);
           scores(i) = mean(diag(r(perm, :)));       
       end
       [best_score, best] = max(scores);
       best_perm  = all_perms(best,:);
   else
       best_perm = (1:num_types) * Hungarian(-r) ;
       best_perm = [best_perm, setdiff(1:num_types, best_perm)];
       best_score = mean(diag(r(best_perm, :)));  
   end

   % fprintf('Best mean corr is %g\n', best_score);
  
   H = H(best_perm, :);
   W = W(:, best_perm);
   
   r_proportions = corr(W, GT_proportions','type',corr_type);
   proportions_score = mean(diag(r_proportions));
     
   
end


function comb=allcomb(ip)
        ncells=length(ip);
        [nd{1:ncells}]=ndgrid(ip{:});
        catted=cat(ncells,nd{1:ncells});
        comb=reshape(catted,length(catted(:))/ncells,ncells);
end
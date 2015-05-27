function [W, H, best_score, proportions_score, ind_scores] = ...
        match_profiles_to_gt(W, H, GT, GT_proportions, corr_type, ...
                             check_for_best_perm)
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
   if ~exist('check_for_best_perm','var')
       check_for_best_perm = true;
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
       ind_scores = best_score;
       best_perm  = all_perms(best,:);
   else
       if num_types_GT==1
           [best_score, best_perm] = max(r);
           ind_scores = best_score;
       else
           best_perm = (1:num_types) * Hungarian(-r) ;
           best_perm = [best_perm, setdiff(1:num_types, best_perm)];

    %        best_score = mean(diag(r(best_perm, :)));  
           best_score = mean_corr_coeff(diag(r(best_perm, :)));  
           ind_scores = diag(r(best_perm, :));
       end
   end

   % fprintf('Best mean corr is %g\n', best_score);
  
   H = H(best_perm, :);
   W = W(:, best_perm);
   
   
%    r_proportions = corr(W, GT_proportions','type',corr_type);
%    proportions_score = mean(diag(r_proportions));

   num_types_W = min(num_types_GT, size(W,2));
   r_proportions = KLDiv(W(:,1:num_types_W), GT_proportions');
   proportions_score = mean(r_proportions);
   
end


function comb=allcomb(ip)
        ncells=length(ip);
        [nd{1:ncells}]=ndgrid(ip{:});
        catted=cat(ncells,nd{1:ncells});
        comb=reshape(catted,length(catted(:))/ncells,ncells);
end
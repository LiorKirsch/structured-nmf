function [W, H, best_score] = match_profiles_to_gt(W, H, GT)
%
% Find a permutation of the rows of H that matches best the
% ground-truth data GT
%
   num_types = size(H, 1);  
   num_types_GT = size(GT, 1); 
   
   if (num_types  < num_types_GT )
       all_perms = allcomb(repmat({1:num_types},num_types_GT,1));
   else
       all_perms = perms(1:num_types);
   end
   num_perms = size(all_perms,1);   
   
   r = corr(H', GT');
   
   scores = zeros(num_perms,1);
   for i=1:num_perms
       perm = all_perms(i,:);
       scores(i) = mean(diag(r(perm, :)));       
   end
   [best_score, best] = max(scores);
   best_perm  = all_perms(best,:);
   % fprintf('Best mean corr is %g\n', best_score);
   
   H = H(best_perm, :);
   W = W(:, best_perm);
   
end


function comb=allcomb(ip)
        ncells=length(ip);
        [nd{1:ncells}]=ndgrid(ip{:});
        catted=cat(ncells,nd{1:ncells});
        comb=reshape(catted,length(catted(:))/ncells,ncells);
end
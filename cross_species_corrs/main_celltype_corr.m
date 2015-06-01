parms.a=1;

load('mouse_cell_type_profiles.mat');
addpath('/cortex/code/cellmix/evaluation/');
addpath('/cortex/code/cellmix/visualization/');




grp =unique(reference);
grp(3) = [] ; % remove Chung which has no cortex samples
calibrated_corr = nan(length(grp));
median_corr = nan(length(grp));
for i =1:length(grp)
    rel_sampl_i = ismember(reference,grp(i)) & is_cortex_or_hippocampus  &is_neuron ;
    rel_sampl_i = logical(sample2type * double(rel_sampl_i));
    
    for j = i+1: length(grp)
        rel_sampl_j = ismember(reference,grp(j)) & is_cortex_or_hippocampus;
        rel_sampl_j = logical(sample2type * double(rel_sampl_j));
        
        corr_score = corr(expression(:,rel_sampl_i), ...
                          expression(:,rel_sampl_j),'type','spearman');
        calibrated_corr(j,i) = mean_corr_coeff(corr_score(:));
        median_corr(j,i) = median(corr_score(:));
    end
end
% 
% figure;
% imagescwithnan(calibrated_corr,hot,[.94 .94 .94]) %# [0 1 1] is cyan
%   
%   
%     ax = gca;
%     ax.YTick = 1:length(grp);
%     ax.YTickLabel = grp;
%     ax.XTick = 1:length(grp);
%     ax.XTickLabel = grp;
%     ax.XTickLabelRotation	=45;
    
    
figure;
imagescnan(median_corr,'NanColor',[.94 .94 .94])
colormap(hot); colorbar;

    ax = gca;
    ax.YTick = 1:length(grp);
    ax.YTickLabel = grp;
    ax.XTick = 1:length(grp);
    ax.XTickLabel = grp;
    ax.XTickLabelRotation	=45;
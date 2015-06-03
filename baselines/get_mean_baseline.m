function [baseline_celltype_profile, baseline_proportions] = ...
        get_mean_baseline(curr_X,gene_inds_predictions)
%
     baseline_celltype_profile = cellfun(@(x) mean(x), curr_X,'UniformOutput',false);
     baseline_celltype_profile = cellfun(@(x) x(:,gene_inds_predictions)', ...
                                         baseline_celltype_profile,'UniformOutput',false);
     baseline_celltype_profile = cellfun(@(x) repmat(x,1,3), ...
                                         baseline_celltype_profile,'UniformOutput',false);
     baseline_proportions = cellfun(@(x) zeros(3,3), curr_X ,'UniformOutput',false);
end
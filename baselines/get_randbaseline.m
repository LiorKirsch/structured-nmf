function result = get_randbaseline(curr_X, gene_inds_predictions,region_names,...
    true_profiles,K,parms)
% Pick random samples from the data set as the profiles 
% and use these samples to match the ground-true profiles.
%

    num_repeats = 300;
    rng(parms.rand_seed);
    
    randbase_score = nan(num_repeats,1);
    randbase_celltype_region_avg_scores = nan(num_repeats,K);
    randbase_region_scores = nan(num_repeats,length(region_names));
%     randbase_celltype_score = nan(num_repeats,1);
    
    for i=1:num_repeats
            printPercentCounter(i, num_repeats);

          [randbase_celltype_profile, randbase_proportions] =...
           get_randsamp_baseline(curr_X,gene_inds_predictions,parms);
     randbase_result = get_all_scores(randbase_celltype_profile, randbase_proportions, ...
                                        true_profiles, region_names,false);  
                                    
                                    
         randbase_score(i) = randbase_result.run_score;
         randbase_celltype_region_avg_scores(i,:) = randbase_result.celltype_region_avg_scores';
         randbase_region_scores(i,:) = randbase_result.region_scores' ;
%          randbase_celltype_score(i) = randbase_result.celltype_scores ;
         
    end

    result.run_score = mean(randbase_score);
    result.run_score_sem = std(randbase_score) /sqrt(num_repeats);
    result.celltype_region_avg_scores = mean(randbase_celltype_region_avg_scores,1);
    result.celltype_region_avg_scores_sem = std(randbase_celltype_region_avg_scores,1,1) /sqrt(num_repeats);
    result.region_scores = mean(randbase_region_scores,1); 
    result.region_scores_sem = std(randbase_region_scores,1,1) /sqrt(num_repeats); 
    result.celltype_scores  = nan;
    result.regions = region_names;
    
end

function [baseline_celltype_profile, baseline_proportions] =...
           get_randsamp_baseline(curr_X,gene_inds_filter,parms)
     
     [num_samp, num_genes] = cellfun(@size, curr_X);
     rand_inds = arrayfun(@(x) randperm(x, parms.num_types), num_samp,'UniformOutput',false);
     
     baseline_celltype_profile = cellfun(@(x,samp_inds) x(samp_inds,gene_inds_filter)',...
         curr_X,rand_inds,'UniformOutput',false);
%      baseline_celltype_profile = cellfun(@(x) x(:,gene_inds_predictions)', ...
%                                          baseline_celltype_profile,'UniformOutput',false);
     baseline_proportions = cellfun(@(x) zeros(3,3), curr_X ,'UniformOutput',false);
end

function inds = choosing_k_inds(n,k, num_repeats)

    all_comb = nchoosek(n,k);
    inds = nan(num_repeats,k);
    
    if all_comb < num_repeats
        inds = nchoosek(1:n,k);
    else
        for i=1:num_repeats
            inds(i,:) = randperm(n,k);
        end
    end
end

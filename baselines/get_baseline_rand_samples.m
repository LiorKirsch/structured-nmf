function [mean_score, std_score, mean_proportions_score, std_proportions_score] = ...
    get_baseline_rand_samples(expression, true_profiles,true_proportions,num_profiles,parms)
% Pick random samples from the data set as the profiles 
% and use these samples to match the ground-true profiles.
%

    num_tissues = size(expression,1);
    num_repeats = 300;

    mean_score = nan(length(num_profiles),1);
    std_score = nan(length(num_profiles),1);
    mean_proportions_score = nan(length(num_profiles),1);
    std_proportions_score = nan(length(num_profiles),1);
    
    for j = 1:length(num_profiles)
        curr_num_profile = num_profiles(j);
        randscores = nan(num_repeats,1);
        rand_proportions_score = nan(num_repeats,1);
        for i=1:num_repeats
            samples = ceil(rand(1,curr_num_profile)*num_tissues);
            HH = expression(samples,:);
            WW = get_proportion_from_profile(expression,HH, parms);
            [~, ~, randscores(i),rand_proportions_score(i)] = match_profiles_to_gt(...
                WW, HH, true_profiles, true_proportions,parms.corr_type);
        end

        mean_score(j) = mean(randscores);
        std_score(j) = std(randscores);
        mean_proportions_score(j) = mean(rand_proportions_score);
        std_proportions_score(j) = std(rand_proportions_score);
    end
end

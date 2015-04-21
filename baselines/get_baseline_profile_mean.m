function [mean_score, std_score, mean_proportions_score, std_proportions_score] = get_baseline_profile_mean(expression, true_profiles,true_proportions,num_profiles,parms)
% Compute the mean profile of all given samples.
% create k new profile which are the mean(expression) + rand*2*std(expression)
% and use these samples to match the ground-true profiles.
%

    mean_profile = mean(expression,1);
    profile_std = std(expression,1,1);
    
    num_tissues = size(expression,1);
    num_repeats = 300;

    mean_score = nan(length(num_profiles),1);
    std_score = nan(length(num_profiles),1);
    mean_proportions_score = nan(length(num_profiles),1);
    std_proportions_score = nan(length(num_profiles),1);

    for j = 1:length(num_profiles)
        curr_num_profile = num_profiles(j);
        
        curr_mean_profile = repmat(mean_profile,[curr_num_profile,1]);
        curr_profile_std = repmat(profile_std,[curr_num_profile,1]);
        
        randscores = nan(num_repeats,1);
        rand_proportions_score = nan(num_repeats,1);
        for i=1:num_repeats
            HH = curr_mean_profile + 2*curr_profile_std.* (2*rand(size(curr_profile_std)) -1);
            WW = get_proportion_from_profile(expression,HH, parms);
            [~, ~, randscores(i),rand_proportions_score(i)] = match_profiles_to_gt(WW, HH, true_profiles, true_proportions,parms.corr_type);
        end

        mean_score(j) = mean(randscores);
        std_score(j) = std(randscores);
        mean_proportions_score(j) = mean(rand_proportions_score);
        std_proportions_score(j) = std(rand_proportions_score);
    end
    
end

function [score, proportions_score] = get_baseline_priors(expression, true_profiles,true_proportions,num_profiles,prior_profiles,parms)
% Use the prior profile as the target profile.
% With it, extract the proportions and compute the correlations.
%
   
    HH = prior_profiles;
    WW = get_proportion_from_profile(expression,HH, parms);
    [~, ~, score,proportions_score] = match_profiles_to_gt(WW, HH, true_profiles, true_proportions,parms.corr_type);
    
end

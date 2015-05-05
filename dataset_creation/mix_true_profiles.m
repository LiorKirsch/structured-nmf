function [proportions,expression, base_proportions] = ...
    mix_true_profiles(true_profiles, region, parms)
% Using different proportion for each region (taken from a survey), we
% create a mix of the true profiles.
% We also multiply each sample with:  randn * (measurement noise)
%
num_tissues = parms.num_tissues;
proportion_variability = parms.proportion_variability;
measurement_noise = parms.measurement_noise;

switch region
    % http://journal.frontiersin.org/article/10.3389/fnana.2014.00127/full
    case {'cortex', 'cortex_L5A', 'cortex_L5B', 'cortex_L6'}
        p_neuro = 0.7; 
        p_astro = 0.1;
        p_oligo = 0.2;
    case 'cerebellum'
        p_neuro = 0.5; 
        p_astro = 0.15;
        p_oligo = 0.35;
    case {'striatum', 'brainstem', 'spinal_cord'}
        p_neuro = 0.65; 
        p_astro = 0.1;
        p_oligo = 0.25;
    otherwise
        error('unknown region');
end
base_proportions =  [p_neuro, p_astro, p_oligo];

expression = zeros(size(true_profiles,1), num_tissues);
proportions = zeros(num_tissues, 3);

for i_tissue = 1:num_tissues
    if mod(i_tissue,100)==0, fprintf('.');end
    p = base_proportions .* (randn(1,3)*proportion_variability+1);
    p = p/sum(p);
    mix = p * true_profiles';
    expression(:,i_tissue) = mix .* max(0,(randn(size(mix))*measurement_noise+1));
    proportions(i_tissue, 1:3) = p;
end

fprintf('\n');
end

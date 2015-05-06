
%
% main script for training/testing cell mix deconvolution
%

% TODO
%
% 3. repeat #2 for other than Cahoy
% 4. What to do when all of the values in the proportion vector are zeros
%       -   add priors to the proportions
%       -   add some small noise to push it a bit from zero - rand(n,1) *eps;
%       -   

init;
parms = conf(parms);

mix_dir = '/cortex/data/microarray/mouse/Okaty2011/Mixtures/';
mix_files = {
    'okaty2011-doyle_cortex_l5a_MN0.01_PR65-10-25_PVAR0.1.mat',...
    'okaty2011-doyle_cortex_l5b_MN0.01_PR65-10-25_PVAR0.1.mat',...
    'okaty2011-doyle_cortex_l6_MN0.01_PR65-10-25_PVAR0.1.mat',...
    'okaty2011-doyle_striatum_MN0.01_PR65-10-25_PVAR0.1.mat',...
    'okaty2011-doyle_cerebellum_MN0.01_PR50-15-35_PVAR0.1.mat',...
    'okaty2011-doyle_brainstem_MN0.01_PR65-10-25_PVAR0.1.mat',...
    'okaty2011-doyle_spinal_cord_MN0.01_PR65-10-25_PVAR0.1.mat'};
parms.mix_files = mix_files;
mix_data = create_multi_region_mix(fullfile(mix_dir, mix_files));
% draw_dendogram(mix_data.profiles, mix_data.cell_types, mix_data.region,'pearson');

parms.structre_type = 'tree';
[parms.structure_matrix,parms.tree_regions] = get_tree_structure(true,mix_data.region);

% parms.structre_type = 'relations';


alg_list = {'alsPinv', 'alsActiveSet', 'mm'}; % 'alsBlockpivot','cjlin', 'prob'}; 
num_samples_list = 10;% [5, 10, 20, 50, 100,200];
num_type_list = 3 ;%1:8;
W_constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};
H_lambda_list = [inf 0 0.001 0.01 0.1 1 10 100 1000];

loop_over_var_name = {};
loop_over_var_value = {};
loop_over_var_name{end + 1} = 'W_constraints'; 
loop_over_var_value{end + 1} = W_constraints_list;
loop_over_var_name{end + 1} = 'num_samples';
loop_over_var_value{end + 1} = num_samples_list;
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;

num_samples =  cellfun( @(x) size(x,1), mix_data.expression);
rand_subset = createRandSubset(num_samples, num_samples_list, parms.subsample_repeats);


X = mix_data.expression;
GT_profiles = mix_data.profiles;
GT_proportions = mix_data.proportions;
parms.regions = mix_data.region;
parms.cell_types = mix_data.cell_types;

if parms.log_transform
    GT_profiles = log2(GT_profiles);
    X = log2(X);
end

%rergion_relatedness_matrix1 = create_r...m(sample_group, profile_group)
%rergion_relatedness_matrix12 = create_r...m(sample_group, profile_group)
% loop_over_var_name{end + 1} = 'RM';
% loop_over_var_value{end + 1} = rm_list;

% [baseline_scores, baseline_std] = get_baseline_profile_mean(X, GT_profiles,loop_over_var_value{end});
[scores, proportions_scores] = loopOverHyperparms_multi(X, GT_profiles, ...
    GT_proportions, rand_subset, parms, ...
    loop_over_var_name, loop_over_var_value, '');
%%

new_loop_over_var_name = loop_over_var_name(2:end);
new_loop_over_var_value = loop_over_var_value(2:end);
for i = 1:length(loop_over_var_value{1});
    
    fprintf('figure - %s\n', loop_over_var_value{1}{i});
    parms.(loop_over_var_name{1}) = loop_over_var_value{1}{i};
    cur_scores = scores{i};
    cur_proportions_scores = proportions_scores{i};
    
    parms.fig_x_axis = new_loop_over_var_name{end};
    figure('Name',sprintf('scores for %s' , loop_over_var_value{1}{i}));
    

    plot_with_var(new_loop_over_var_name, new_loop_over_var_value, cur_scores,cur_proportions_scores,X,GT_profiles, GT_proportions, parms);
    subplot(1,2,1);
    title(sprintf('profile - %s',loop_over_var_value{1}{i}));
    subplot(1,2,2);
    title(sprintf('proportions - %s',loop_over_var_value{1}{i}));
    ylim([-1 1 ]);
    
    fig_file_name = set_filenames('figure', parms);
    saveas(gcf,fig_file_name);
   
end
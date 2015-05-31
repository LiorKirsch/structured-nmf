
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
mix_files = {};
mix_files = [ mix_files ,    {'okaty2011-doyle_cortex_l5a_MN0.2_PR70-10-20_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_cortex_l5b_MN0.2_PR70-10-20_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_cortex_l6_MN0.2_PR70-10-20_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_striatum_MN0.2_PR65-10-25_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_cerebellum_MN0.2_PR50-15-35_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_brainstem_MN0.2_PR65-10-25_PVAR0.2.mat'}];
mix_files = [ mix_files ,    {'okaty2011-doyle_spinal_cord_MN0.2_PR65-10-25_PVAR0.2.mat'}];

parms.mix_files = mix_files;
mix_data = create_multi_region_mix(fullfile(mix_dir, mix_files));
% draw_dendogram(mix_data.profiles, mix_data.cell_types, mix_data.region,'pearson');
% 
% parms.structre_type = 'tree';


% parms.structre_type = 'relations_dist';
% parms.structre_type = 'relations_parentdist';
parms.structre_type = 'relations_parent_level';
% parms.structre_type = 'relations_on_expression';

[parms.structure_matrix,parms.relation_regions] = get_tree_structure(mix_data.region);
[parms.structure_matrix,parms.relation_regions] = get_relation_structure(parms.structure_matrix,parms.relation_regions,mix_data.region,parms.structre_type,mix_data.expression);

parms.do_sep_init = true;
parms.num_samples = 5; 

% alg_list = {'alsPinv', 'alsActiveSet', 'mm'}; % 'alsBlockpivot','cjlin', 'prob'}; 
alg_list = {'alsActiveSet'}; 
num_samples_list = [ 5, 10, 20, 50];%, 100,200];
num_type_list = 3 ;%1:8;
W_constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};
W_constraints_list = { 'on_simplex_with_noise'};

% W_constraints_list = {'on_simplex_with_noise'};
% H_lambda_list = [  1 10 ];
H_lambda_list = [ 0 0.001 0.01 0.1 1 10 100 1000 inf];
% H_lambda_list = 1;
H_lambda_list = [100 1000 inf];

parms.num_restarts = 30; 
parms.subsample_repeats = 30; 
parms.init_type = 'random';


loop_over_var_name = {};
loop_over_var_value = {};
loop_over_var_name{end + 1} = 'W_constraints'; 
loop_over_var_value{end + 1} = W_constraints_list;
% loop_over_var_name{end + 1} = 'num_markers';
% loop_over_var_value{end + 1} = [5,20,50,100];
loop_over_var_name{end + 1} = 'num_samples';   % Must appear as part of loop
loop_over_var_value{end + 1} = num_samples_list; % Must appear as part of loop
loop_over_var_name{end + 1} = 'nmf_method';
loop_over_var_value{end + 1} = alg_list;
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
% [scores,proportions_scores, parmslist] = iterateParms_multi(X, GT_profiles, ...
%     GT_proportions, rand_subset, parms, ...
%     loop_over_var_name, loop_over_var_value, true);

[scores, proportions_scores] = loopOverHyperparms_multi(X, GT_profiles, ...
    GT_proportions, rand_subset, parms, ...
    loop_over_var_name, loop_over_var_value, '');
%%
parms.draw_log_scale = true;
new_loop_over_var_name = loop_over_var_name(2:end);
new_loop_over_var_value = loop_over_var_value(2:end);


colors_to_use = linspecer(length(loop_over_var_value{2}));
colors_to_use = autumn(length(loop_over_var_value{2}));
set(groot,'defaultAxesColorOrder',colors_to_use);


for i = 1:length(loop_over_var_value{1});

    
    if iscell(loop_over_var_value{1})
        curr_val_string = loop_over_var_value{1}{i};
        parms.(loop_over_var_name{1}) = loop_over_var_value{1}{i};
    else
        curr_val_string = sprintf('%g', loop_over_var_value{1}(i));
        parms.(loop_over_var_name{1}) = loop_over_var_value{1}(i);
    end
    
    fprintf('figure - %s\n', curr_val_string);
    figure('Name',sprintf('scores for %s' , curr_val_string));
        
    cur_scores = scores{i};
    cur_proportions_scores = proportions_scores{i};
    
    parms.fig_x_axis = new_loop_over_var_name{end};
    
    

    plot_with_var(new_loop_over_var_name, new_loop_over_var_value, cur_scores,cur_proportions_scores,X,GT_profiles, GT_proportions, parms);
    ylim([.88 0.95]);
    subplot(1,2,1);
    title(sprintf('profile - %s',curr_val_string));
    subplot(1,2,2);
    title(sprintf('proportions - %s',curr_val_string));
    ylim([0 1.5 ]);
    
    fig_file_name = set_filenames('figure', parms);
    saveas(gcf,fig_file_name);
   
end
set(groot,'defaultAxesColorOrder','default');       
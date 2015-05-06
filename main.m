
%
% main script for training/testing cell mix deconvolution
%

% TODO
% 2. Add regularizer weights... 
%
% 3. repeat #2 for other than Cahoy
% 4. What to do when all of the values in the proportion vector are zeros
%       -   add priors to the proportions
%       -   add some small noise to push it a bit from zero - rand(n,1) *eps;
%       -   

init;

nmf_method = take_from_struct(parms, 'nmf_method', 'alsBlockpivot');


% Load the mixture data
parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.1_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.05_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.01_PR60-10-30_PVAR0.1';

% parms.dataset_file = 'barres2014-lin-lin_MN0.1_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'barres2014-lin-lin_MN0.01_PR60-10-30_PVAR0.1';

mix_data = load_data(parms.dataset_file, parms);

% Load the priors
% parms.prior_dataset = 'Doyle';
% parms.prior_types = {'neuro', 'astro', 'oligo'};
parms.prior_dataset = 'Cahoy';
parms.prior_types = {'neuro', 'astro', 'oligo'};

prior_data = load_priors(parms.prior_dataset, parms.prior_types, parms);

parms.num_types = 3;
parms.num_samples = 50;
parms.maxiter = 500;
parms.loglevel = 0;
parms.W_constraints = 'on_simplex';
% parms.corr_type = 'Spearman';
parms.corr_type = 'Pearson';

parms.record_scores = true;
parms.rand_seed = 42; % The answer to life the universe and everything
parms.num_restarts = 5; % <===  increase to 30

parms.H_lambda = 0.1;
parms.H_prior = prior_data; 

parms.W_lambda = 0;

parms.log_transform = false;

parms.subsample_repeats = 5; % <=== increase to 30 
alg_list = {'alsPinv', 'alsActiveSet', 'mm'}; % 'alsBlockpivot','cjlin', 'prob'}; 
alg_list = {'alsActiveSet'}; % 'alsBlockpivot','cjlin', 'prob'}; 
num_samples_list = 10;% [5, 10, 20, 50, 100,200];
num_type_list = 3 ;%1:8;
W_constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};
% W_constraints_list = {'on_simplex'};
% W_constraints_list = {'inside_simplex'};
W_constraints_list = {'positive'};
H_lambda_list = [0 0.001 0.01 0.1 1 10 100 1000];
loop_over_var_name = {};
loop_over_var_value = {};

loop_over_var_name{end + 1} = 'W_constraints';
loop_over_var_value{end + 1} = W_constraints_list;
loop_over_var_name{end + 1} = 'nmf_method';
loop_over_var_value{end + 1} = alg_list;
loop_over_var_name{end + 1} = 'num_samples';
loop_over_var_value{end + 1} =  num_samples_list;
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;





[num_genes,num_samp] = size(mix_data.expression);
rand_subset = createRandSubset(num_samp, num_samples_list, parms.subsample_repeats);

X = mix_data.expression';
GT_profiles = mix_data.profiles';
GT_proportions = mix_data.proportions';

if parms.log_transform
    GT_profiles = log2(GT_profiles);
    X = log2(X);
end

% [baseline_scores, baseline_std] = get_baseline_profile_mean(X, GT_profiles,loop_over_var_value{end});
[scores, proportions_scores] = loopOverHyperparms(X, GT_profiles, GT_proportions, rand_subset,parms, ...
                            loop_over_var_name, loop_over_var_value, '');
%%

new_loop_over_var_name = loop_over_var_name(2:end);
new_loop_over_var_value = loop_over_var_value(2:end);
for i = 1:length(loop_over_var_value{1} );
    
    
    fprintf('figure - %s\n', loop_over_var_value{1}{i} );
    parms.(loop_over_var_name{1}) = loop_over_var_value{1}{i};
    cur_scores = scores{i};
    cur_proportions_scores = proportions_scores{i};
    
    parms.fig_x_axis = new_loop_over_var_name{end};
    figure('Name',sprintf('scores for %s' , loop_over_var_value{1}{i} ));
    

    plot_with_var(new_loop_over_var_name, new_loop_over_var_value, cur_scores,cur_proportions_scores,X,GT_profiles, GT_proportions, parms);
    subplot(1,2,1);
    title(sprintf('profile - %s',loop_over_var_value{1}{i} ));
    subplot(1,2,2);
    title(sprintf('proportions - %s',loop_over_var_value{1}{i} ));
    ylim([-1 1 ]);
    
    fig_file_name = set_filenames('figure', parms);
    saveas(gcf,fig_file_name);
   
end

%
% main script for training/testing cell mix deconvolution
%

% TODO:
% 1. Clean and save into GIT
% 2. Figure 1: 
%    mean-corr as a function number of samples. 
%    repeat for the various methods
%    (draw several sampels of size 5)
%
% 3. Add regularizer weights... 
%
%
% 4. repeat #2 for other than Cahoy, also other noise levels etc.
%
%
%

init;

nmf_method = take_from_struct(parms, 'nmf_method', 'alsBlockpivot');


% Load the mixture data
parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.01_PR60-10-30_PVAR0.1';
mix_data = load_data(parms.dataset_file, parms);


parms.num_types = 3;
parms.num_samples = 50;
parms.maxiter = 500;
parms.loglevel = 0;
% parms.W_constrains = 'on_simplex';
% parms.W_constrains = 'inside_simplex';
parms.W_constrains = 'positive';

parms.record_scores = true;
parms.rand_seed = 42; % The answer to life the universe and everything
parms.num_restarts = 10;

parms.use_regularizer = false;
parms.H_lambda = 100;
parms.W_lambda = 0;

parms.subsample_repeats = 10;
alg_list = {'alsPinv','alsBlockpivot','alsActiveSet','mm'}; % 'cjlin', 'prob'}; 
num_samples_list =50 ;% [5, 10, 20, 50, 100,200];
num_type_list = 1:8;
loop_over_var_name = {};
loop_over_var_value = {};

loop_over_var_name{end + 1} = 'nmf_method';
loop_over_var_value{end + 1} = alg_list;
loop_over_var_name{end + 1} = 'num_samples';
loop_over_var_value{end + 1} =  num_samples_list;
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;

[num_genes,num_samp] = size(mix_data.expression);
rand_subset = createRandSubset(num_samp, num_samples_list, parms.subsample_repeats);

X = mix_data.expression';
GT_profiles = mix_data.profiles';
scores = loopOverHyperparms(X,GT_profiles,rand_subset,parms, loop_over_var_name, loop_over_var_value ,'');


parms.fig_x_axis = loop_over_var_name{end};
plot_with_var(loop_over_var_name, loop_over_var_value, scores);
fig_file_name = set_filenames('figure', parms);
saveas(gcf,fig_file_name);
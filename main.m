
%
% main script for training/testing cell mix deconvolution
%

% TODO:
% 1. Figure 1: 
%    mean-corr as a function number of samples. 
%    repeat for the various methods
%    (draw several sampels of size 5)
%
% 2. Add regularizer weights... 
%
%
% 3. repeat #2 for other than Cahoy, also other noise levels etc.
%
%
%

init;

nmf_method = take_from_struct(parms, 'nmf_method', 'alsBlockpivot');


% Load the mixture data
% parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.1_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.05_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-lin_cahoy_MN0.01_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-log_cahoy_MN0.1_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'okaty2011-lin-log_cahoy_MN0.05_PR60-10-30_PVAR0.1';

parms.dataset_file = 'barres2014-lin-lin_MN0.1_PR60-10-30_PVAR0.1';
% parms.dataset_file = 'barres2014-lin-lin_MN0.01_PR60-10-30_PVAR0.1';

mix_data = load_data(parms.dataset_file, parms);


parms.num_types = 3;
parms.num_samples = 50;
parms.maxiter = 500;
parms.loglevel = 0;
parms.W_constrains = 'on_simplex';

parms.record_scores = true;
parms.rand_seed = 42; % The answer to life the universe and everything
parms.num_restarts = 10;

parms.H_lambda = 0;
parms.W_lambda = 0;

parms.subsample_repeats = 10;
alg_list = {'alsPinv','alsActiveSet','mm'}; % 'alsBlockpivot','cjlin', 'prob'}; 
num_samples_list =50; %[5, 10, 20, 50, 100,200];
num_type_list = 8:9;
% W_constrains_list = {'on_simplex', 'inside_simplex', 'positive'};
W_constrains_list = {'on_simplex'};
% W_constrains_list = {'inside_simplex'};
% W_constrains_list = {'positive'};
loop_over_var_name = {};
loop_over_var_value = {};

loop_over_var_name{end + 1} = 'W_constrains';
loop_over_var_value{end + 1} = W_constrains_list;
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
%%

new_loop_over_var_name = loop_over_var_name(2:end);
new_loop_over_var_value = loop_over_var_value(2:end);
for i = 1:length(loop_over_var_value{1} );
    fprintf('figure - %s\n', loop_over_var_value{1}{i} );
    parms.(loop_over_var_name{1}) = loop_over_var_value{1}{i};
    cur_scores = scores{i};
    
    parms.fig_x_axis = new_loop_over_var_name{end};
    figure('Name',sprintf('scores for %s' , loop_over_var_value{1}{i} ));
    plot_with_var(new_loop_over_var_name, new_loop_over_var_value, cur_scores);
    fig_file_name = set_filenames('figure', parms);
    saveas(gcf,fig_file_name);
end
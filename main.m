
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
mix_data = load_data('okaty2011_cahoy_MN0.1_PR60-10-30_PVAR0.1', parms);


parms.num_types = 4;
parms.num_samples = 50;
parms.maxiter = 500;
parms.loglevel = 0;
% parms.W_constrains = 'on_simplex';
parms.W_constrains = 'inside_simplex';
% parms.W_constrains = 'positive';

parms.record_scores = true;
parms.rand_seed = 42; % The answer to life the universe and everything
parms.num_restarts = 10 ;

parms.use_regularizer = false;
parms.H_lambda = 100;
parms.W_lambda = 0;

alg_list = {'alsPinv','alsBlockpivot','alsActiveSet','mm'}; % 'cjlin', 'prob'}; 
% alg_list = {'alsPinv'}; 
% alg_list = {'alsBlockpivot'}; 
% alg_list = {'alsActiveSet'};
% alg_list = {'mm'}; 
[num_genes, num_samp] = size(mix_data.expression);

num_samples_list = [5, 10, 20, 50, 100,200];
num_types_list = 1:8;
subsample_repeats = 10;

% rand_subset = createRandSubset(num_samp, num_samples_list, subsample_repeats);
rand_subset = createRandSubset(num_samp, parms.num_samples, subsample_repeats);

alg_scores = {};
for i_alg = 1:length(alg_list)
   parms.nmf_method = alg_list{i_alg};
   
%     scores = nan(length(num_samples_list),subsample_repeats);
%     for i_ns = 1:length(num_samples_list)
%         fprintf('===== samples size %d =====\n', num_samples_list(i_ns) );
%         parms.num_samples = num_samples_list(i_ns);
%         curr_samples_selected = rand_subset{i_ns};
        
    scores = nan(length(num_types_list),subsample_repeats);
    for i_ns = 1:length(num_types_list)
        fprintf('===== num types %d =====\n', num_types_list(i_ns) );
        parms.num_types = num_types_list(i_ns);
        curr_samples_selected = rand_subset{1};

        for j_sr = 1:subsample_repeats
            current_parms = parms; % to activiate parfor
            current_parms.subsample_iter = j_sr;
            samples_selected = curr_samples_selected(j_sr,:);
            X = mix_data.expression(:,samples_selected)';
            [W, H, diff_record, time_record, eucl_dist] = load_nmf_results(X, parms.num_types, current_parms.nmf_method, current_parms);

            % Match profiles to ground truth
            [W, H, best_score] = match_profiles_to_gt(W,H, mix_data.profiles');
            fprintf('Best mean corr is %g\n', best_score);    

            scores(i_ns,j_sr) = best_score;
        end
        % Evaluate
    end
    alg_scores{i_alg} = scores;
end


% plot_with_var(num_samples_list, alg_scores, alg_list);
plot_with_var(num_types_list, alg_scores, alg_list);
  
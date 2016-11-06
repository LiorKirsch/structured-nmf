
%
% main script for comparing different methods of cellmix deconvolution
%


init;
parms = conf(parms);

mix_dir = '/cortex/data/microarray/mouse/Okaty2011/Mixtures/';
mix_files = {};
mix_files = [ mix_files ,    {'okaty2011-doyle_cortex_l5a_MN0.2_PR70-10-20_PVAR0.2.mat'}];


parms.mix_files = mix_files;
parms.mix_dir = mix_dir;
parms = rmfield(parms, 'structre_type');
parms.do_sep_init = true;
parms.num_samples = 5; 
parms.num_restarts = 5; 
parms.subsample_repeats = 5; 
parms.init_type = 'random';
parms.init_subtype = 'random';
parms.maxiter = 3000;


alg_list = {'alsAccProj','alsBlockpivot','alsActiveSet','accProj','mm','cjlin'}; % 'alsPinv' 'alsBlockpivot','cjlin', 'prob'}; 
% alg_list = {'alsActiveSet'}; 
num_samples_list = [  5,10 20, 50];%, 100,200];
num_type_list = 3 ; % 1:8; % [3,4,5] ;
W_constraints_list = {'on_simplex', 'inside_simplex', 'positive'} ;%,'on_simplex_with_noise'};
% W_constraints_list = {'on_simplex'} ;
H_lambda_list = 0;


loop_over_var_name = {};
loop_over_var_value = {};
loop_over_var_name{end + 1} = 'W_constraints'; 
loop_over_var_value{end + 1} = W_constraints_list;
% loop_over_var_name{end + 1} = 'num_markers';
% loop_over_var_value{end + 1} = [5,20,50,100];
loop_over_var_name{end + 1} = 'nmf_method';
loop_over_var_value{end + 1} = alg_list;
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;
loop_over_var_name{end + 1} = 'num_samples';   % Must appear as part of loop
loop_over_var_value{end + 1} = num_samples_list; % Must appear as part of loop
                

% num_samples =  size(mix_data.expression,1);
% rand_subset = createRandSubset(num_samples, num_samples_list, parms.subsample_repeats);
% rand_subset = rand_subset{1};

is_synthetic = true;
[results, H_results, W_results] = loopOverHyperparms_new(parms, loop_over_var_name, ...
                                     loop_over_var_value , 'rand_run',is_synthetic);

%%
parms.draw_log_scale = false;
new_loop_over_var_name = loop_over_var_name(2:end);
new_loop_over_var_value = loop_over_var_value(2:end);



for i = 1:length(loop_over_var_value{1});

    set(0,'defaultAxesFontName', 'Dejavu Sans')
    set(0,'defaultTextFontName', 'Dejavu Sans')
    if iscell(loop_over_var_value{1})
        curr_val_string = loop_over_var_value{1}{i};
        parms.(loop_over_var_name{1}) = loop_over_var_value{1}{i};
    else
        curr_val_string = sprintf('%g', loop_over_var_value{1}(i));
        parms.(loop_over_var_name{1}) = loop_over_var_value{1}(i);
    end
    
    fprintf('figure - %s\n', curr_val_string);
    FigHandle = figure('Name',sprintf('scores for %s' , curr_val_string));
    set(gca,'FontSize', 18);
    
    cur_scores = scores{i};
    cur_proportions_scores = proportions_scores{i};
    
    parms.fig_x_axis = new_loop_over_var_name{end};
    
    plot_with_var(new_loop_over_var_name, new_loop_over_var_value, cur_scores,cur_proportions_scores,X,GT_profiles, GT_proportions, parms);

    fig_file_name = set_filenames('figure', parms);
    saveas(gcf,fig_file_name);
    print(gcf,'-depsc ',[fig_file_name(1:end-4),'.eps']);
end

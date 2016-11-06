% load dermantitis cell types

% % ====== GPL15520 array (138 samples, 22088 genes, 3 adults, 3 prenatal)
% load('/cortex/data/RNA-Seq/human/Darmanis2015/rnaseq_celltypes_GPL15520.mat');

% ====== GPL18573 array (328 samples, 22088 genes, 5 adults, 1 prenatal)
% load('/cortex/data/RNA-Seq/human/Darmanis2015/rnaseq_celltypes_GPL18573.mat');



init;
parms = conf(parms);

[~, parms] = take_from_struct(parms, 'dataset', 'brainspan2014');
% parms.dataset = 'human6';
[default_regions, parms.species] = get_region_set(parms.dataset);
[~, parms] = take_from_struct(parms, 'regions_ack', default_regions);
parms.regions_ack = sort(parms.regions_ack);

parms.true_dataset = 'darmanis';


% parms = rmfield(parms, 'structre_type');
parms.do_sep_init = true;
parms.num_restarts = 5; 
parms.init_type = 'random';
parms.init_subtype = 'random';
parms.maxiter = 3000;


% alg_list = {'alsAccProj','alsBlockpivot','alsActiveSet','accProj','mm','cjlin'}; % 'alsPinv' 'alsBlockpivot','cjlin', 'prob'}; 
alg_list = {'alsAccProj'}; 
num_samples_list = 5 ; %[  5,10 20, 50];%, 100,200];
num_type_list = 3 ; % 1:8; % [3,4,5] ;
% W_constraints_list = {'on_simplex', 'inside_simplex', 'positive'} ;%,'on_simplex_with_noise'};
W_constraints_list = {'on_simplex'} ;
H_lambda_list = 0;
H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0, 10.^[-3:1:3], inf] );

loop_over_var_name = {};
loop_over_var_value = {};
loop_over_var_name{end + 1} = 'H_lambda';
loop_over_var_value{end + 1} = H_lambda_list;
loop_over_var_name{end + 1} = 'num_samples';   % Must appear as part of loop
loop_over_var_value{end + 1} = num_samples_list; % Must appear as part of loop
loop_over_var_name{end + 1} = 'W_constraints'; 
loop_over_var_value{end + 1} = W_constraints_list;
% loop_over_var_name{end + 1} = 'num_markers';
% loop_over_var_value{end + 1} = [5,20,50,100];
loop_over_var_name{end + 1} = 'nmf_method';
loop_over_var_value{end + 1} = alg_list;
loop_over_var_name{end + 1} = 'num_types';
loop_over_var_value{end + 1} = num_type_list;

                             
[results, H_results, W_results, parms_results] = loopOverHyperparms_new(parms, loop_over_var_name, ...
    loop_over_var_value , 'collect_results');



draw_celltype_scores(loop_over_var_name, loop_over_var_value, ...
                           results, parms, 'Corr');

% if take_from_struct(parms, 'do_plot', true);
%     parms.draw_log_scale = true;
%     % draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
% 
%     switch getenv('USER')
%         case 'lior', 
%           draw_indv_figure(loop_over_var_name, loop_over_var_value, ...
%                            results, parms, 'Corr');
%       case 'gal', 
%           draw_indv_figure_gal(loop_over_var_name, loop_over_var_value, ...
%                                results, parms, 'Corr');
%     end
% end

[best_score, base_score] = report_results(results, parms);
    

 
 
 
 
 
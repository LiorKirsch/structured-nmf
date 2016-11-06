init
parms = conf(parms);
parms.do_plot = false;
parms.do_save = false;
parms.regions_short = {'M1C', 'V1C'};
parms.regions_short = {'DFC', 'V1C'};
parms.regions_short = {'DFC', 'HIP'};
parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'V1C'};

parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'STC' 'V1C'};

parms.regions_short =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
                    ...
                    'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};


% parms.regions_short = {'DFC' 'OFC', 'MFC', 'VFC'};
parms.init_subtype = 'noise10';

parms.true_dataset = 'barres';
% parms.W_constraints = 'positive';
parms.num_restarts = 50;
parms.max_iter = 1000;


parms.H_lambda_list=[0 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 ...
                    1000];

init
parms = conf(parms);
parms.do_plot = false; 
parms.do_save = false; 
parms.regions_short = {'M1C', 'V1C'};
parms.regions_short = {'DFC', 'V1C'};
parms.regions_short = {'DFC', 'HIP'};
parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'V1C'};

parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'STC' 'V1C'};

parms.regions_short =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
                    'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};


% parms.regions_short = {'DFC' 'OFC', 'MFC', 'VFC'};
parms.init_subtype = 'noise10';

parms.true_dataset = 'barres';
% parms.W_constraints = 'positive';
parms.num_restarts = 50;
parms.max_iter = 1000;


parms.H_lambda_list=[0 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000]; 

% parms.H_lambda_list=[0 0.1 1 10 100]; 
main

% draw_neuron_figure(parms.H_lambda_list, results, parms);

draw_indv_figure_gal(loop_over_var_name, loop_over_var_value, results, parms);


region_string = sprintf('%s-', parms.regions_short{1:end});
fig_filename = sprintf('Figures/neuron_%s', region_string);

region_string = sprintf('%s-', parms.regions_short{1:end});
region_string = region_string(1:end-1);
fig_filename = sprintf('Figures/3panels_var_%s', region_string);


do_save = take_from_struct(parms, 'do_save', false);  
fig_save(fig_filename, do_save);


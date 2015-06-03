init
parms = conf(parms);
parms.do_plot=false; 
parms.regions_short = {'M1C', 'V1C'};
parms.regions_short = {'DFC', 'V1C'};
parms.regions_short = {'DFC', 'HIP'};
parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'V1C'};
parms.init_subtype='noise10';

parms.H_lambda_list=[0 0.1 1 5 10 50 100]; 
main

draw_neuron_figure(parms.H_lambda_list, results, parms);
region_string = sprintf('%s-', parms.regions_short{1:end});
fig_filename = sprintf('Figures/neuron_%s', region_string);

do_save = take_from_struct(parms, 'do_save', false);  
fig_save(fig_filename, do_save);


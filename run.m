init
parms = conf(parms);
parms.do_plot = false;
parms.init_type = 'random';
parms.init_subtype = 'noise10';


parms.regions_short = {'M1C', 'S1C'};
parms.regions_short = {'CBC', 'ITC'};
parms.do_plot=1;
parms.true_dataset = 'barres';


parms.H_lambda_list = [0 0.01 0.1 0.1 1 10]; main



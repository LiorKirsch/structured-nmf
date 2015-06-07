init
parms = conf(parms);
parms.do_plot = false; 
parms.do_save = false; 
parms.regions_short = {'M1C', 'V1C'};
parms.regions_short = {'DFC', 'V1C'};
parms.regions_short = {'DFC', 'HIP'};
parms.regions_short = {'AMY' 'CBC', 'DFC', 'HIP', 'V1C'};

parms.init_subtype = 'noise10';
parms.true_dataset = 'barres';
parms.num_restarts = 50;
parms.max_iter = 1000;

parms.H_lambda_list=1; 
main

H_results

CALL = cell2mat(H_results{1}{1}');


CDFC = CALL(7:9,:)';
CV1C = CALL(13:15,:)';
CCBC = CALL(4:6,:)';

G = CCBC;

inds = find((G(:,1)+1) ./ (G(:,2)+G(:,3)+1)>5)'
G(inds,:)
gene_info.gene_symbols{inds}

%CN = CALL(1:3:15,:)';
%inds = find(CN(:,2) - CN(:,3)>4)
%gene_info.gene_symbols{inds}
%[CN(inds,2) - CN(inds,3)]









region_string = sprintf('%s-', parms.regions_short{1:end});
fig_filename = sprintf('Figures/neuron_%s', region_string);

region_string = sprintf('%s-', parms.regions_short{1:end});
region_string = region_string(1:end-1);
fig_filename = sprintf('Figures/3panels_var_%s', region_string);


do_save = take_from_struct(parms, 'do_save', false);  
fig_save(fig_filename, do_save);


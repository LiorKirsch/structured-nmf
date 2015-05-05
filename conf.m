function parms = conf(parms)


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
parms.W_lambda = 0;

parms.log_transform = false;

parms.subsample_repeats = 5; % <=== increase to 30 


end
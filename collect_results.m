

parms.dataset = 'brainspan2014'; 
      
all_regions =  {'A1C','CBC','DFC','HIP', 'IPC','ITC', ...
             'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'}; %'AMY'

% regions =  {'DFC','M1C','S1C','V1C'};
% regions =  {'DFC','M1C','V1C'};
% regions =  {'M1C','V1C'};
% regions =  {'M1C','S1C'}; 

% task = 'collect';
task = 'generate';

all_pairs = nchoosek(1:length(all_regions), 2);


pair_best_scores = nan(size(all_pairs,1),1);
pair_base_scores = nan(size(all_pairs,1),1);

for  i = 1:size(all_pairs,1)
      
      parms.regions_short = all_regions(all_pairs(i,:) );
      parms.init_type = 'random';
      parms.init_subtype = 'random';
%       parms.W_constraints = 'positive';
      parms.W_constraints = 'on_simplex_with_noise';
      
%       parms.num_restarts = 30;
      parms.num_restarts = 5;
%       parms.maxiter = 1000;  
      parms.maxiter = 500;  
%       parms.H_lambda_list = [0, 10.^[-3:0.5:3], inf];
      parms.H_lambda_list = [ 0 0.001 0.01 0.1 1 10 100 1000 inf];


    switch task
        case 'collect'
            parms.only_collect = true;
            parms.do_plot = false;
              main
              disp(all_regions(all_pairs(i,:)) );
              pair_best_scores(i) = best_score;
              pair_base_scores(i) = base_score;
              fprintf('=============== file found for');
        case 'generate'
            parms.only_collect = false;
            parms.do_plot = false;
            main
    end

end


[best_diff,best_ind] = nanmax(pair_best_scores - pair_base_scores);
fprintf('best score %g for :', best_diff);
disp( all_regions(all_pairs(best_ind,:) ) );


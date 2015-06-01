

parms.dataset = 'brainspan2014'; 
      
region_set =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
             'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};

rng('shuffle') ;
all_pairs = nchoosek(1:length(region_set), 2);
all_pairs = all_pairs(randperm(size(all_pairs,1)),:);
for  i = 1:size(all_pairs,1)
%     try
%     parms.structure_filter = 10;
      
      parms.regions = region_set(all_pairs(i,:) );
      parms.init_type = 'random';
      parms.W_constraints = 'positive';
      parms.num_restarts = 30;
      parms.maxiter = 1000;  
      parms.H_lambda_list = [0, 10.^[-3:0.5:3], inf];

      
      disp('=========================================');
      disp(parms.regions);
      disp('=========================================');
      
       main_real_data 
%     catch
        
%     end
end

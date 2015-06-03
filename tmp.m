

parms.dataset = 'brainspan2014'; 
      
regions =  {'A1C','AMY','CBC','DFC','HIP', 'IPC','ITC', ...
             'M1C','MFC','OFC','S1C', 'STC','V1C','VFC'};

% regions =  {'DFC','M1C','S1C','V1C'};
% regions =  {'DFC','M1C','V1C'};
% regions =  {'M1C','V1C'};
% regions =  {'M1C','S1C'}; 


rng('shuffle') ;
all_pairs = nchoosek(1:length(regions), 2);
all_pairs = all_pairs(randperm(size(all_pairs,1)),:);
for  i = 1:size(all_pairs,1)
%     try
%     parms.structure_filter = 10;
      
      parms.regions_short = regions(all_pairs(i,:) );
      parms.init_type = 'random';
      parms.W_constraints = 'positive';
      parms.init_subtype =  'samples_with_noise';
      
      parms.num_restarts = 50;
      parms.maxiter = 1000;  
      parms.H_lambda_list = [0, 10.^[-3:0.5:3], inf];

      
      disp('=========================================');
      disp(parms.regions_short);
      disp('=========================================');
      
       main
%     catch
        
%     end
end

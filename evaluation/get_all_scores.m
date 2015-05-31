
function output = get_all_scores(predicted_profiles, predicted_proportions, ...
                                               true_profiles, ...
                                               region_name, show_individual)
   if ~exist('show_individual','var')
       show_individual = true;
   end
   
    num_regions = length(predicted_profiles);
    assert(length(true_profiles) == num_regions, ...
           'true and predicted profiles should have the same number of elements');        
    best_scores = nan(num_regions,1);
    individual_scores = cell(num_regions,1);
    
    
    for i = 1:num_regions
        % need to transpose the expression
        GT_proportions = zeros(size(predicted_proportions{i},1), ...
                               size(true_profiles{i},2));
        [~, ~, best_scores(i), ~,individual_scores{i}] = ...
            match_profiles_to_gt(predicted_proportions{i}, ...
                                 predicted_profiles{i}', true_profiles{i}', ...
                                 GT_proportions', 'spearman'); 
        
        % when the number of true cell type is less then 3 pad the vector with nans
        individual_scores{i} = [individual_scores{i}; ...
                nan( max(0,3 -length(individual_scores{i})),1)]; 
            
%         fprintf('%25s: %4.2f ',region_name{i}, 100*best_scores(i));
%         
%         if show_individual
%             for individual_i = 1:length(individual_scores{i})
%                 fprintf('\ttype %d: %4.2f', individual_i, 100*individual_scores{i}(individual_i));
%             end
%             fprintf('\n');
%         end
    end
    output.celltype_scores =  individual_scores;
    output.region_scores =  best_scores;
    output.regions =  region_name;        
    output.run_score = mean_corr_coeff(best_scores);

    dim = ndims(individual_scores);          % Get the number of dimensions for your arrays
    M = cat(dim+1,individual_scores{:});        % Convert to a (dim+1)-dimensional matrix
    output.celltype_region_avg_scores = nanmean(M,dim+1);     % Get the mean across arrays
                
%     fprintf('\tREGION MEAN SCORE: %5.3g====\n', output.run_score);
end
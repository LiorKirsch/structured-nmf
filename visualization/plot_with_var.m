function plot_with_var(loop_over_var_name, loop_over_var_value, scores)
% scores as a matrix of size [num_sizes x repeats]
% sample_sizes as a vector of size [num_sizes x 1]

    figure('Name','corr vs sample_size');
    hold on;

    alg_name = loop_over_var_value{1};
    
    assert( strcmp('nmf_method', loop_over_var_name{1} ),'nmf_method should be the first');
    
    x = loop_over_var_value{end};
    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    for i = 1:length(alg_name)
        %changing scores to - [different_parm X repeats]
        curr_scores = scores{i};
        for j=2: (length(loop_over_var_value) -1)
           curr_scores = curr_scores{1}; 
        end

        % turn the last value to the last cell
        curr_scores = cell2mat(curr_scores);
        curr_scores = curr_scores';
        
        
        y = mean(curr_scores,2);
        y_std = std(curr_scores,1,2);
        errorbar(x,y,y_std,'LineWidth',2);
    end
    legend(alg_name,'Location','southeast');
    hold off;
%     xlabel('num types');
    xlabel(x_label);
    ylabel('correlation');
    ylim([0.85 1]);
end
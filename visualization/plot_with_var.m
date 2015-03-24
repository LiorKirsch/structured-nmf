function plot_with_var(sample_sizes, scores, alg_name)
% scores as a matrix of size [num_sizes x repeats]
% sample_sizes as a vector of size [num_sizes x 1]

    figure('Name','corr vs sample_size');
    hold on;
    x = sample_sizes;
    for i = 1:length(alg_name)
        y = mean(scores{i},2);
        y_std = std(scores{i},1,2);
        errorbar(x,y,y_std,'LineWidth',2);
    end
    legend(alg_name,'Location','southeast');
    hold off;
    xlabel('num types');
    ylabel('correlation');
    ylim([0.4 1.1]);
end
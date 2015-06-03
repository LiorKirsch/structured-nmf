function draw_neuron_figure(H_lambda_list, results, parms)
%
    % init
    clrs = set_colors();
    LW = 'Linewidth', CLR = 'Color'; FS = 'fontsize';

    [real_scores, rand_scores] = results_struct_to_mat(results);
    % scores has dims: [num_lambdas, num_regions, num_types]
    xxx = H_lambda_list;
    if numel(xxx) ~= size(real_scores, 1)
        error();
    end
    xxx(xxx == 0) = 10^-5;
    xxx(isinf(xxx)) = 10^+5;

    scores = squeeze(real_scores(:,:,1)); % neuron data only
    rands = squeeze(rand_scores(:,:,1));
    figure(1); clf; hold on;


    y_type = take_from_struct(parms, 'y_type', 'variance');    
    switch y_type
      case 'corr',     ylabel = 'spearman';
      case 'variance', ylabel = 'unexplained variance';
        scores = 1- scores.^2;
    end
    
    y_mmm = mean(scores, 2);
    y_sss = std(scores, [], 2);
    r_mmm = mean(rands, 2);
    r_sss = std(rands, [], 2);        

    
    errorbar(xxx, y_mmm, y_sss, '-', CLR, clrs{1}, LW, 3);
    plot(x, y_mmm, '-', CLR, clrs{1}, LW, 3);
    set(gca, 'xscale', 'log')
    
    set(gca, 'xtick', x);
    xlabel('\lambda', FS, 20)
    ylabel(ylabel, FS, 20)

    
end
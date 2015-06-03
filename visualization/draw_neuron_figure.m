function draw_neuron_figure(H_lambda_list, results, parms)
%
    % init
    clrs = set_colors();
    LW = 'Linewidth'; CLR = 'Color'; FS = 'fontsize';

    [real_scores, rand_scores] = results_struct_to_mat(results);
    % scores has dims: [num_lambdas, num_regions, num_types]
    xxx = H_lambda_list;
    if numel(xxx) ~= size(real_scores, 1)
        error();
    end
    
    lambda0_ind = find(xxx==0);
    if isempty(lambda0_ind) 
        error('\n\n\t\tERROR: missing HL=%d\n\n\n', 0); 
    end
    iii = find(xxx>0);
    xxx = xxx(iii);

    scores = squeeze(real_scores(:,:,1)); % neuron data only
    rands = squeeze(rand_scores(:,:,1));
    figure(1); clf; hold on;

    % Transform the values if needed
    y_type = take_from_struct(parms, 'y_type', 'variance');    
    switch y_type
      case 'corr',     ylbl = 'spearman';
      case 'variance', ylbl = 'unexplained variance';
        scores = 1- scores.^2;
        rands = 1- rands.^2;        
    end
    y_mmm = mean(scores, 2)'
    y_sss = std(scores, [], 2)' / sqrt(2)
    r_mmm = mean(rands, 2)'
    r_sss = std(rands, [], 2)' / sqrt(2);
    b_mmm = y_mmm(lambda0_ind)*ones(size(xxx))
    b_sss = y_sss(lambda0_ind)*ones(size(xxx)) / sqrt(2);
    
    % Plot results
    h1 = errorbar(xxx, y_mmm(iii), y_sss(iii), '-', CLR, clrs{1}, LW, 2);
    h2 = errorbar(xxx, r_mmm(iii), r_sss(iii), '-', CLR, clrs{3}, LW, 2);
    h3 = errorbar(xxx, b_mmm, b_sss, '-', CLR, clrs{2}, LW, 2);
    set(gca, 'xscale', 'log');
    set(gca, 'ytick', 0.30:0.01:0.6);    

    % Title 
    region_string = sprintf('%s, ', parms.regions_short{1:end});
    title(sprintf('regions = %s', region_string(1:end-2)));
    
    xlabel('\lambda', FS, 20);
    ylabel(ylbl, FS, 20);
end
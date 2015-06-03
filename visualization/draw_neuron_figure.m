function draw_neuron_figure(loop_over_var_name, loop_over_var_value, ...
                            results, parms, y_label)
%

    [real_scores, rand_scores] = results_struct_to_mat(results);
    % scores has dims: [num_lambdas, num_regions, num_types]
    x = loop_over_var_value{2};
    if numel(x) ~= size(real_scores, 1)
        error();
    end
    
    if take_from_struct(parms, 'draw_log_scale', true)
       x(x == 0) = 10^-5;
       x(isinf(x)) = 10^+5;
    end
    
    % For plotting neuron data, we take the first cell type only
    scores = squeeze(real_scores(:,:,1));
    rands = squeeze(rand_scores(:,:,1));
    
    clrs = set_colors();
    y = mean(scores, 2); y = y(:)';
    figure(1); clf; hold on;
    plot(x, y, '-', 'Color', clrs{0})
    
    
    
    
    curr_results = results{1};
      cell_type=1;
XXX      for j=1:3
            assert(length(curr_results) ==  length(x),'problem with build the results');
            y = cellfun(@(x) x.celltype_region_avg_scores(j), curr_results);
            subplot(1,3,j);hold on;
            plot_h(i) = plot(x,y);
        end
    end
    
    subplot(1,3,1); title('neuro');
    subplot(1,3,2); title('astro');
    subplot(1,3,3); title('oligo');
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%s %g', ...
                                               loop_over_var_name{1},x), ...
                                  legend_strings, 'UniformOutput', false);
    end
    

    set(groot,'defaultAxesColorOrder','default');       
    % draw the mean baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(1,3,j);hold on;
            baseline = cellfun(@(x) x.baseline_celltype_region_avg_scores(j), curr_results);
            plot_h(n +1) = plot(x,baseline,'--');
        end
        
    end
    legend_strings{end +1} = sprintf('mean-profile baseline');
    

    set(groot,'defaultAxesColorOrder','default');       
    % draw the rand baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(1,3,j);hold on;
            baseline = cellfun(@(x) ...
                               x.randbase_celltype_region_avg_scores(j), ...
                               curr_results);
            baseline_sem = cellfun(@(x) ...
                                   x.randbase_celltype_region_avg_scores_sem(j), ...
                                   curr_results);
            plot_h(n +2) = errorbar(x,baseline,baseline_sem,'.');
        end
    end
    legend_strings{end +1} = sprintf('rand-sample baseline');
    legend_strings = strrep(legend_strings, '_', ' ');
    legend(plot_h, legend_strings, 'location', 'best');
    x_label = strrep(loop_over_var_name{end}, '_', ' ') ;
    xlabel(x_label);
    ylabel(y_label);
    
    if parms.draw_log_scale
        for j=1:3
            subplot(1,3,j); hold on;
            set(gca,'XScale','log');
            xticks = 10.^[-3:1:3];
            set(gca,'XTick', xticks);
            tmp = get(gca, 'XTickLabel');
            if any(loop_over_var_value{end} == 0 )
                tmp{loop_over_var_value{end} == 0} = '0';
            end
            if any(isinf(loop_over_var_value{end}))
                tmp{isinf(loop_over_var_value{end})} = 'inf';
            end
            set(gca,'XTickLabel', log10(xticks));
        end
    end
    % suptitle(strrep(set_parmstr(parms), '_', ' '));
    header = sprintf('%s-', parms.regions_short{1:end});
    header = sprintf('%s %s:%s', header(1:end-1), parms.init_type, parms.init_subtype);
    header = sprintf('%s NN%d', header, parms.num_restarts);
    suptitle(header);
    parms.fig_x_axis = x_label;
    [~,file_name,dir_name] = set_filenames('figure_real', parms);
    fprintf('drawing figure %s\n',fullfile(dir_name,file_name));
    saveas(gcf,fullfile(dir_name,['indv_', file_name]));
end



function check_argins(loop_over_var_value, results);    
    recursion_levels = length(loop_over_var_value);
    msg = 'results should have 2 levels';
    assert(recursion_levels == 2 , msg);
    legend_strings = loop_over_var_value{1};
    msg = 'problem with build the results';
    assert(length(results) == length(legend_strings), msg);
end
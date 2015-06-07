function draw_indv_figure_gal(loop_over_var_name, loop_over_var_value, ...
                              results, parms, y_label)
%
    set(gca,'FontSize', 18);
    num_rows = 1; 

    switch getenv('USER')
      case 'lior', 
        set(FigHandle, 'Position', [100, 100, 1600, 1000]);
        FigHandle = figure('Name',parms.dataset_file); clf; hold on;
      case 'gal', 
        figure(1); clf; hold on;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition',[0 0 10 15]);
    end
    set(groot, 'defaultAxesColorOrder', 'default');
    
    n = length(loop_over_var_value{1});    
    check_argins(loop_over_var_value, results);    
    legend_strings = loop_over_var_value{1};
    x = loop_over_var_value{2};    

    if take_from_struct(parms, 'draw_log_scale')
       x(x == 0) = 10^-4;
       x(isinf(x)) = 10^+4;
    end

    plot_h = nan(n,1);
    for i = 1:n
        curr_results = results{i};
        for j=1:3
            assert(length(curr_results) ==  length(x),'problem with build the results');
            y = cellfun(@(x) x.celltype_region_avg_scores(j), curr_results);
            subplot(num_rows,3,j);hold on;
            y_var = (1-y.^2);
            plot_h(i) = plot(x, y_var, 's-');
        end
    end    
    subplot(num_rows,3,1); title('neuron');
    subplot(num_rows,3,2); title('astrocyte');
    subplot(num_rows,3,3); title('oligodendrocyte');
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%s %g', ...
                                               loop_over_var_name{1},x), ...
                                  legend_strings, 'UniformOutput', false);
    end

    set(groot,'defaultAxesColorOrder','default');       
    % draw the rand baselines
    for i = 1:n
        curr_results = results{i};
         for j=1:3
            subplot(num_rows,3,j);hold on;
            set(gca,'XScale','log')
            baseline = cellfun(@(x) ...
                               x.randbase_celltype_region_avg_scores(j), ...
                               curr_results);
            baseline_sem = cellfun(@(x) ...
                                   x.randbase_celltype_region_avg_scores_sem(j), ...
                                   curr_results);
            baseline_var = 1-baseline.^2;
            plot_h(n +1) = errorbar(x, baseline_var, baseline_sem*4 ,'--');

        end
    end
    
    do_legend = take_from_struct(parms, 'do_legend', false)
    if do_legend
        legend_strings{end +1} = sprintf('matched samples');
        legend_strings = strrep(legend_strings, '_', ' ');    
        legend_strings{1} = 'inferred profiles';    
        legend(plot_h, legend_strings, 'location', 'best', 'box', 'off');
    end

    subplot(num_rows,3,1);
    ylabel('spearman correaltion');
    ylabel('unexplained variance');
    
    for j=1:3
        hs(j) = subplot(num_rows,3,j); hold on;        
        xlabel('log_{10} \lambda');
        set(gca,'xtick', [10^-4 10^-2 1 10^2 10^4]);
        set(gca,'xticklabels', {'ind' '-2' '0' '2' '4' });


        if 0
            switch x_label
              case '\lambda'
                xlim([10^-5.5, 10^5.5]);
                if parms.draw_log_scale 
                    skip_array = [1,2,4,6,8];
                    x_ticks = x(skip_array);
                    set(gca,'XTick', x_ticks);
                    tmp = get(gca,'XTickLabel');
                    tmp = strrep(tmp, '10^{-4}','indv') ;
                    tmp = strrep(tmp, '10^{4}','joined') ;
                    set(gca,'XTickLabel', tmp);
                end
            end
        end


    end

    do_save = take_from_struct(parms, 'do_save', true)
    if do_save
        parms.fig_x_axis = strrep(loop_over_var_name{end}, ' ', '_');    
        [~, file_name, dir_name] = set_filenames('figure_real', parms);
        fprintf('save figure [%s]\n', fullfile(dir_name,file_name));
        saveas(gcf,fullfile(dir_name,['indv_', file_name]));
    end
end



function check_argins(loop_over_var_value, results);    
    recursion_levels = length(loop_over_var_value);
    msg = 'results should have 2 levels';
    assert(recursion_levels == 2 , msg);
    legend_strings = loop_over_var_value{1};
    msg = 'problem with build the results';
    assert(length(results) == length(legend_strings), msg);
end
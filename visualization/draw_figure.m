function draw_figure(loop_over_var_name, loop_over_var_value, results, parms,y_label)
%
    figure('Name',parms.dataset_file); clf; hold on;
    recursion_levels = length(loop_over_var_value);
    
    assert(recursion_levels ==2 , ...
           'results should have 2 levels - the first for the legend, the second for the x axis');
    
    legend_strings = loop_over_var_value{1};
    x = loop_over_var_value{2};
    
    assert(length(results) == length(legend_strings), 'problem with build the results');

    if parms.draw_log_scale
       x(x == 0) = 10^-7;
       x(isinf(x)) = 10^+7;
    end

    
    for i = 1:length(legend_strings)
        curr_results = results{i};
        assert(length(curr_results) ==  length(x),'problem with build the results');
        y = cellfun(@(x) x.run_score, curr_results);
        plot(x,y);
    end
    
    
    if ~iscell(legend_strings)
        legend_strings = arrayfun(@(x) sprintf('%g',x), legend_strings, ...
                                  'UniformOutput', false);
    end
    

    % draw the mean baselines
    for i = 1:length(legend_strings)
        baseline = cellfun(@(x) x.baseline_score, curr_results);
        plot(x,baseline,'--');
        legend_strings{end +1} = sprintf('mean profile (%s)', legend_strings{i});
    end

    
        
    legend_strings = strrep(legend_strings, '_', ' ');
    legend(legend_strings, 'location', 'best');
    x_label = strrep(loop_over_var_name{end},'_',' ') ;
    xlabel(x_label);
    ylabel(y_label);
    
    if parms.draw_log_scale
        set(gca,'XScale','log')
        set(gca,'XTick', x);
        tmp = get(gca,'XTickLabel');
        if any(loop_over_var_value{end} == 0 )
            tmp{loop_over_var_value{end} == 0} = '0';
        end
        if any(isinf(loop_over_var_value{end}))
            tmp{isinf(loop_over_var_value{end})} = 'inf';
        end
        set(gca,'XTickLabel', tmp);
    end
    
        

end
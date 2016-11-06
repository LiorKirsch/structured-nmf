function draw_celltype_scores(loop_over_var_name, loop_over_var_value, results, parms,y_label)
% The first variable in loop_over_... is used for the x-axis
% The scoend variable is used for the legends

    
    
    set(groot,'defaultAxesColorOrder','default');       
    
    x = loop_over_var_value{1};    

    if take_from_struct(parms, 'draw_log_scale')
       x(x == 0) = 10^-5;
       x(isinf(x)) = 10^+5;
    end

    % plot the scores:
    
    scores = unpack_scores_from_struct('celltype_region_avg_scores', ...
        loop_over_var_name, loop_over_var_value, results, parms);
    baseline = unpack_scores_from_struct('baseline_celltype_region_avg_scores', ...
        loop_over_var_name, loop_over_var_value, results, parms);
%     baseline_sem = unpack_scores_from_struct('randbase_celltype_region_avg_scores_sem', ...
%         loop_over_var_name, loop_over_var_value, results, parms);
    randbase = unpack_scores_from_struct('randbase_celltype_region_avg_scores', ...
        loop_over_var_name, loop_over_var_value, results, parms);
    randbase_sem = unpack_scores_from_struct('randbase_celltype_region_avg_scores_sem', ...
        loop_over_var_name, loop_over_var_value, results, parms);
    
    
    % loop over the cell types neuron,astro,oligo
    cell_types = {'neuro','astro','oligo'};
    for k = 1:3 
        cell_scores = cellfun(@(x) x(k), scores);
        cell_baseline = cellfun(@(x) x(k), baseline);
%         cell_baseline_sem = cellfun(@(x) x(k), baseline_sem);
        cell_randbase = cellfun(@(x) x(k), randbase);
        cell_randbase_sem = cellfun(@(x) x(k), randbase_sem);
        subplot(1,3,k);
        plot(x, cell_scores, 's-');
        hold on;
        plot(x, cell_baseline);
        errorbar(x,cell_randbase,cell_randbase_sem,'--');
        title(cell_types{k});
    end
    
    legend_strings = {'reconstructed matched','samples-mean baseline','random baseline'};
    legend(legend_strings, 'location', 'best');
    
    subplot(1,3,1);
    ylabel(y_label);
    
    for k=1:3
        subplot(1,3,k);
        set(gca,'XScale','log')
        xlim([10^-5.5, 10^5.5]);
        if parms.draw_log_scale 
%             skip_array = [1,2,4,6,8,9];
%             x_ticks = x(skip_array);
            x_ticks = x;
            set(gca,'XTick', x_ticks);
            tmp = get(gca,'XTickLabel');
            tmp = strrep(tmp, '10^{-5}','indv') ;
            tmp = strrep(tmp, '10^{5}','joined') ;
            set(gca,'XTickLabel', tmp);
        end
    end
%     parms.fig_x_axis = strrep(loop_over_var_name{end}, ' ', '_');
%     [~, file_name, dir_name] = set_filenames('figure_real', parms);
%     fprintf('save figure [%s]\n', fullfile(dir_name,file_name));
%     saveas(gcf,fullfile(dir_name,['indv_', file_name]));
end



function check_argins(loop_over_var_value, results)    
    recursion_levels = length(loop_over_var_value);
    msg = 'results should have 2 levels';
    assert(recursion_levels == 2 , msg);
    legend_strings = loop_over_var_value{1};
    msg = 'problem with build the results';
    assert(length(results) == length(legend_strings), msg);
end
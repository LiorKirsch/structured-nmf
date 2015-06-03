init;
parms = conf(parms);
only_collect = take_from_struct(parms, 'only_collect', false);

if only_collect
    X = nan;
    gene_info = nan;
else

    dataset = take_from_struct(parms, 'dataset', 'brainspan2014');
    [default_regions, parms.species] = get_region_set(dataset);
    regions_short = take_from_struct(parms, 'regions_short', default_regions);
    regions_short = sort(regions_short);
    
    % Load all data
    [parms.dataset_file, expression, gross_region_vec, gene_info, ...
     gross_structures_info] = load_all(dataset, regions_short);
    
    % limit to a set of specific genes
    expression = change_to_linear_scale(expression);
    
    gross_regions = gross_structures_info(gross_region_vec);
    [X, region_names] = split_to_cell(expression, gross_regions);
    clear('expression', 'gross_region_vec', 'gross_regions', ...
          'gross_structures_info');
    
    fprintf('main: Get the structure matrix\n');
    parms.regions = region_names;
    parms = get_structure_matrix(parms.dataset_file, parms.structre_type, ...
                                               region_names, parms);
end

    num_type_list = take_from_struct(parms, 'num_type_list', [3]) ;%1:8;
    
% %     Removed markers from file name
% %     num_markers_list = take_from_struct(parms, 'num_markers_list', [5 20 50 100 200]);
    %%%
    
    H_lambda_list = take_from_struct(parms, 'H_lambda_list', [0, 10.^[-3:0.1:3], inf] );
    gene_subset_list = take_from_struct(parms, 'gene_subset_list', ...
                                        {'all','okaty_anova10000', ...
                        'okaty_anova5000' 'okaty_infogain5000', ...
                        'okaty_infogain10000', 'okaty_infogain1000'});
    gene_okaty_filter_list = take_from_struct(parms, 'gene_okaty_filter_list', ...
                                              {'cortex_doyle', 'cortex', 'doyle', 'all'});
    constraints_list = {'on_simplex', 'inside_simplex', 'positive','on_simplex_with_noise'};
    H_lambda_list = sort(H_lambda_list);

    loop_over_var_name = {};
    loop_over_var_value = {};
    loop_over_var_name{end + 1} = 'num_types';
    loop_over_var_value{end + 1} = num_type_list;
    loop_over_var_name{end + 1} = 'H_lambda';
    loop_over_var_value{end + 1} = H_lambda_list;
    
    parms.cell_types = arrayfun(@(idx) sprintf('#%d',idx),1:parms.num_types,'UniformOutput',false);
    
    fprintf('main: loop over hyper parameters\n');
    results = loopOverHyperparms_real(X, gene_info, parms, ...
                                      loop_over_var_name, ...
                                      loop_over_var_value ,'');
    
    if take_from_struct(parms, 'do_plot', true);
        parms.draw_log_scale = true;
        % draw_figure(loop_over_var_name, loop_over_var_value, results, parms, 'Corr');
        draw_indv_figure(loop_over_var_name, loop_over_var_value, ...
                         results, parms, 'Corr');
    
    end


[best_score, base_score] = report_results(results, parms);
    

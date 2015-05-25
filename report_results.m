function  report_results(results)
%

    num_vars = length(results);
    for i_var = 1: num_vars
        r = results{i_var};        
        num_regions = length(r.regions);

        fprintf('%s = %g\n', r.var_name, r.var_value);
        for i_region = 1:num_regions
            fprintf('\t%20s: ', r.regions{i_region}(1:20));
            fprintf(' %4.2f ', 100*r.celltype_scores{i_region});
            fprintf(' (avg=%4.2f)\n', 100*r.region_scores(i_region));
        end

    end
end

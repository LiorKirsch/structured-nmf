function [best_score, base_score] = report_results(results)

    report_results_struct(results);
    [real_scores, rand_scores] = struct_to_mat(results);
    [num_lambdas, num_regions, num_celltypes] = size(scores);
    [best_score, ind_lambda] = max(squeeze(median(real_scores(:,:,1), 2)));
    base_scores = rand_scores(ind_lambda, 1, 1);
end


function report_results_struct(results)    
%
    num_vars_outer = length(results);
    for i_var_outer = 1: num_vars_outer
        results_outer = results{i_var_outer}; 
        
        num_vars = length(results_outer);
        for i_var = 1: num_vars
            r = results_outer{i_var};        
            num_regions = length(r.regions);

            fprintf('%s  :\n', r.loop_string);
%             fprintf('%s = %g\n', r.var_name, r.var_value);
            for i_region = 1:num_regions
                fprintf('\t%20s: ', r.regions{i_region});
                fprintf(' %4.2f ', 100*r.celltype_scores{i_region});
                fprintf(' (avg=%4.2f)\n', 100*r.region_scores(i_region));
            end
        end
        
        fprintf('==== mean profile baseline ====:\n');
        for i_region = 1:num_regions
            fprintf('\t%20s: ', r.regions{i_region});
            fprintf(' %4.2f ', 100*r.baseline_celltype_score{i_region});
            fprintf(' (avg=%4.2f)\n', 100*r.baseline_region_scores(i_region));
        end
        fprintf('==== rand samples baseline ====:\n');
        for i_region = 1:num_regions
            fprintf('\t%20s: ', r.regions{i_region});
            fprintf(' %4.2f ', 100*r.randbase_celltype_score(i_region));
            fprintf(' (avg=%4.2f)\n', 100*r.randbase_region_scores(i_region));
        end
        fprintf('\t%20s: ', 'mean over regions');
        fprintf(' %4.2f ', 100* mean(r.randbase_celltype_score(i_region)),1);
        fprintf('\n');
        fprintf('===============================:\n');
    end
end

function out = printthis(ob)
    if iscellstr(ob)
        out = sprintf('%s',ob);
    else
        out = sprintf('%g',ob);
    end
end


function [num_lambdas, num_regions, num_celltypes] = get_dims(results)
    num_lambdas = length(results{1});
    r = results{1}{1};
    num_regions = length(r.regions);
    num_celltypes = length(r.celltype_scores{1}); 
end


function [real_scores, rand_scores] = struct_to_mat(results)
    [num_lambdas, num_regions, num_celltypes] = get_dims(results);
    real_scores = zeros(num_lambdas, num_regions, num_celltypes);
    rand_scores = zeros(num_lambdas, num_regions, num_celltypes);
    i_var_outer = 1; 
    results_outer = results{i_var_outer}; 
        
    for i_lambda = 1: num_lambdas
        r = results_outer{i_lambda};        
        for i_region = 1:num_regions
            real_scores(i_lambda, i_region, 1:num_celltypes) = r.celltype_scores{i_region};
            rand_scores(i_lambda, i_region, 1:num_celltypes) = ...
                r.randbase_celltype_scores{i_region};
        end
    end
end
    


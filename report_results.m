function [best_score, base_score] = report_results(results, parms)
%
    report_results_struct(results, parms);
    [real_scores, rand_scores] = results_struct_to_mat(results);
    [num_lambdas, num_regions, num_celltypes] = size(real_scores);
    [best_score, ind_lambda] = max(squeeze(median(real_scores(:,:,1), 2)));
    base_score = rand_scores(ind_lambda, 1, 1);
end


function report_results_struct(results, parms)    
%
    fprintf('\nRESULS FOR');
    fprintf(' %s-', parms.regions_short{1:end});
    fprintf(' (%s:%s)\n', parms.init_type, parms.init_subtype);
    num_vars_outer = length(results);
    for i_var_outer = 1: num_vars_outer
        results_outer = results{i_var_outer}; 
        
        num_vars = length(results_outer);
        for i_var = 1: num_vars
            r = results_outer{i_var};        
            num_regions = length(r.regions);
            header = strrep(r.loop_string, '-', '');
            header = strrep(header, '    ', ' ');
            header = strrep(header, '   ', ' ');
            header = strrep(header, '  ', ' ');
            header = strrep(header, 'num_types 3', '');
            fprintf('    %s\n', header);
            scrs = zeros(num_regions, 3);
            for i_region = 1:num_regions
                region_short = parms.regions_short{i_region};
                region = r.regions{i_region};
                %fprintf('\t%-40s', regions_short):
                %fprintf(' %4.2f ', 100*r.celltype_scores{i_region});
                %fprintf(' (avg=%4.2f)\n', 100* r.region_scores(i_region));
                scrs(i_region,1:3) = r.celltype_scores{i_region};
            end
            fprintf('\t%-40s', 'mean over regions');
            fprintf(' %4.2f ', 100* mean(scrs,1));
            fprintf('\n');
        end
        
        fprintf('     mean profile baseline\n');
        % for i_region = 1:num_regions
        %    fprintf('\t%-40s', r.regions{i_region});
        %    fprintf(' %4.2f ', 100*r.baseline_celltype_score{i_region});
        %    fprintf(' (avg=%4.2f)\n', 100*r.baseline_region_scores(i_region));
        % end
        fprintf('\t%-40s', 'mean over regions');
        fprintf(' %4.2f ', 100* r.baseline_celltype_region_avg_scores);
        fprintf('\n');
        
        fprintf('     rand samples baseline\n');
        % for i_region = 1:num_regions
        %    region_short = parms.regions_short{i_region};
        %    region = r.regions{i_region};        
        %    fprintf('\t%-40s', r.regions{i_region});
        %    fprintf(' %4.2f ', 100*r.randbase_celltype_score{i_region});
        %    fprintf(' (avg=%4.2f)\n', 100*r.randbase_region_scores(i_region));
        % end
        fprintf('\t%-40s', 'mean over regions');
        fprintf(' %4.2f ', 100* r.randbase_celltype_region_avg_scores);
        fprintf('\n');
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



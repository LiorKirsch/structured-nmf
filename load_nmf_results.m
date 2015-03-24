function  [best_W, best_H, best_diff_record, best_time_record, ...
          eucl_dist] = load_nmf_results(X, K, nmf_method, parms)
%
%
%

    vars = {'best_W', 'best_H', 'best_diff_record', ...
            'best_time_record','eucl_dist'};
    filename  = set_filenames('demixing', parms);
    [do_calc, best_W, best_H, best_diff_record, best_time_record, ...
     eucl_dist] = cond_load(filename, 0, vars{1:end});
    if do_calc < 1 
       return
    end
    
    % Compute the demixing
    [best_W, best_H, best_diff_record, best_time_record, eucl_dist] ...
        = nmf(X, K, nmf_method, parms);
        
    save(filename, vars{1:end});
    fprintf('Saved demixing results into [%s]\n', filename);
end
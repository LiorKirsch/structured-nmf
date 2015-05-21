function loopOverHyperparms_real(X, ...
                                     parms, loop_over_var_name, ...
                                     loop_over_var_value ,loop_string)
% This function calls the nmf factorization but with different parms
% it loops over the set of parms provided in 
%  loop_over_var_name
%     and
%  loop_over_var_value
%

if length(loop_over_var_value) ==1
    % The stopping phase of the recursion
    
    var_name = loop_over_var_name{1};
    var_values = loop_over_var_value{1};
    
    parfor i_vars = 1: length(var_values)
       current_parms = parms;
       new_loop_string = add_var_to_loop_string(loop_string, var_values, i_vars, var_name, current_parms);
       set_terminal_title(new_loop_string);
    
        [W, H, diff_record, time_record, eucl_dist] = ...
            load_nmf_results(X, parms.num_types, ...
                             parms.nmf_method, parms);
    end
  
else
    % Remove one layer from the recursion
    var_name = loop_over_var_name{1};
    loop_over_var_name = loop_over_var_name(2:end);
    var_values = loop_over_var_value{1};
    loop_over_var_value = loop_over_var_value(2:end);

    
    for i_vars = 1: length(var_values)
       current_parms = parms;
       new_loop_string = add_var_to_loop_string(loop_string, var_values, i_vars, var_name, current_parms);
       
        loopOverHyperparms_real(X,...
                                current_parms, ...
                                loop_over_var_name, ...
                                loop_over_var_value, ...
                                new_loop_string); 
    end
    
end

set_terminal_title('done');
end


function new_loop_string = add_var_to_loop_string(loop_string, var_values, i_vars, var_name, current_parms)

    if iscellstr(var_values)
       current_parms.(var_name) = var_values{i_vars};
       new_loop_string = sprintf('%s - %s %s',loop_string, var_name, current_parms.(var_name) );
   else
       current_parms.(var_name) = var_values(i_vars);
       new_loop_string = sprintf('%s - %s %g',loop_string, var_name, current_parms.(var_name) );
   end

end
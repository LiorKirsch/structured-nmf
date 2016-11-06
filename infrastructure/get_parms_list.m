function [parms_list,parms_strings] = get_parms_list(parms, loop_over_var_name, ...
                                     loop_over_var_value)
                                 
   [parms_list,parms_strings] = get_parms_list_rec(parms, loop_over_var_name, ...
                                     loop_over_var_value,'',{},{});
                                 
   parms_strings = cellfun(@(x) x(4:end) ,parms_strings,'UniformOutput',false);
end
                                 
function [parms_list,parms_strings] = get_parms_list_rec(parms, loop_over_var_name, ...
                                     loop_over_var_value,loop_string,parms_list,parms_strings)
   if length(loop_over_var_name) ==1
       var_name = loop_over_var_name{1};
       var_values = loop_over_var_value{1};
       new_parms_list = cell(length(var_values),1);
       new_parms_string = cell(length(var_values),1);
       for i_vars =1:length(var_values)
           
            [new_loop_string, current_parms, var_current_value] = ...
             add_var_to_parms(loop_string, var_values, i_vars, ...
                                    var_name, parms);
            
         
          new_parms_list{i_vars} =  current_parms;
          new_parms_string{i_vars} =  new_loop_string;
       end
       parms_list = [new_parms_list; parms_list];
       parms_strings = [new_parms_string; parms_strings];
   else
       
      var_name = loop_over_var_name{1};
      loop_over_var_name = loop_over_var_name(2:end);
      var_values = loop_over_var_value{1};
      loop_over_var_value = loop_over_var_value(2:end);
       
       for i_vars =1:length(var_values)
          [new_loop_string, current_parms, var_current_value] = ...
          add_var_to_parms(loop_string, var_values, i_vars, ...
                                    var_name, parms);
            
          [parms_list,parms_strings] = get_parms_list_rec(current_parms, loop_over_var_name, ...
                                     loop_over_var_value,new_loop_string,parms_list,parms_strings);
          
       end
       
       
   end                           
                                 
                                 
   end

   
   function [new_loop_string, current_parms, value] = ...
        add_var_to_parms(loop_string, var_values, i_vars, ...
                               var_name, current_parms)
    if iscellstr(var_values)
       value = var_values{i_vars};
       current_parms.(var_name) = value;
       new_loop_string = sprintf('%s <> %s %s', loop_string, var_name, ...
                                 current_parms.(var_name) );
    else
       value = var_values(i_vars);
       current_parms.(var_name) = value;
       new_loop_string = sprintf('%s <> %s %5.3g', loop_string, var_name, ...
                                 current_parms.(var_name) );
   end

end
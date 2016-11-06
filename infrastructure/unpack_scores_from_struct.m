function scores = unpack_scores_from_struct(field_name, loop_over_var_name, loop_over_var_value, results, parms)


loop_sizes = cellfun(@length, loop_over_var_value);
loop_var = find(loop_sizes > 1);
assert(length(loop_var) <= 2 & ~isempty(loop_var),'This type of figure takes two running variable');


if length(loop_var) == 1
    
    size1 = loop_sizes(loop_var(1));
    values_1 = loop_over_var_value{loop_var(1)};
    name_1 = loop_over_var_name{loop_var(1)};

    scores = cell(size1 , 1);
    for i=1:length(results)
        [~,ind1] = ismember(results{i}.parms.(name_1), values_1);

        scores{ind1} = results{i}.(field_name);
    end
elseif length(loop_var) == 2
    
    
    size1 = loop_sizes(loop_var(1));
    size2 = loop_sizes(loop_var(2));
    scores = cell(size1 , size2);

    values_1 = loop_over_var_value{loop_var(1)};
    values_2 = loop_over_var_value{loop_var(2)};
    name_1 = loop_over_var_name{loop_var(1)};
    name_2 = loop_over_var_name{loop_var(2)};

    for i=1:length(results)
        [~,ind1] = ismember(parms{i}.(name_1), values_1);
        [~,ind2] = ismember(parms{i}.(name_2), values_2);

        scores{ind1,ind2} = results{i}.(field_name);
    end
end


end
        
    
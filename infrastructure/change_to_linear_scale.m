function expression = change_to_linear_scale(expression)

    if all(expression(:) < 100)
        disp('data seems like it is in log scale - transforming it to linear...');
        expression = 2 .^ expression;
    else
        disp('data seems like it is in linear scale');
    end

end
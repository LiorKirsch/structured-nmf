function expression = change_to_linear_scale(expression)
%
    if all(expression(:) < 100)
        disp('data seems to be in log scale; transforme to linear.');
        expression = 2 .^ expression;
    else
        disp('data seems to be in linear scale.');
    end
end
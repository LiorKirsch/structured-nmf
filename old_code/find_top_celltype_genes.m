
%
% Load the okaty-doyle profiles
%
    init 

    cell_types = load('mouse_cell_type_profiles.mat');
    cell_types.expression = 2.^cell_types.expression;
    cell_type_filter = strmatch('Doyle', cell_types.reference);
    cell_types = limit_data_by_cell_type_filter(cell_types, cell_type_filter);
    exp = cell_types.expression * double(cell_types.sample2type);
    cell_types.expression_types = exp ./ ...
        repmat(sum(cell_types.sample2type,1), size(cell_types.expression,1),1);



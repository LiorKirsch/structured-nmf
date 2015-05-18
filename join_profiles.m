function [joined_profiles, joined_types] = join_profiles(profile_cell, cell_desc)
[num_types, num_features] = size(profile_cell{1});
num_cells = length(profile_cell);

assert(size(cell_desc,1) == num_cells,'each item in the cell should have a description');
joined_profiles = cat(1,profile_cell{:});

joined_types = repmat(cell_desc,1,3)';
joined_types = joined_types(:);

end
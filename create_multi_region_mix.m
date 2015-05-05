function [mix_data] = create_multi_region_mix(mix_files)
% profile_group is a cell array of size:
%      num_cell_types*num_regions x 2
% for example if there are two regions
%     cortex_L5A  pyramidal neuron
%     cortex_L5A  astro
%     cortex_L5A  oligo
%     cortex_L5B  pyramidal neuron
%     cortex_L5B  astro
%     cortex_L5B  oligo

% mix_files = {
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_cortex_l5a_MN0.01_PR65-10-25_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_cortex_l5b_MN0.01_PR65-10-25_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_cortex_l6_MN0.01_PR65-10-25_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_striatum_MN0.01_PR65-10-25_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_brainstem_MN0.01_PR65-10-25_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_cerebellum_MN0.01_PR50-15-35_PVAR0.1.mat',...
%     '/cortex/data/microarray/mouse/Okaty2011/Mixtures/okaty2011-doyle_spinal_cord_MN0.01_PR65-10-25_PVAR0.1.mat'};
%     [profiles,proportions,expression,sample_group,profile_group] = create_multi_region_mix(mix_files);


profiles = cell(length(mix_files),1);
proportions = cell(length(mix_files),1);
expression = cell(length(mix_files),1);
cell_types = cell(length(mix_files),1);
region = cell(length(mix_files),1);
for i = 1:length(mix_files)
    file_name = mix_files{i};
    data = load(file_name);
    
    profiles{i} = data.profiles;
    proportions{i} = data.proportions;
    expression{i} = data.expression';
    cell_types{i} = data.cell_types';
    region{i}  = data.region;
end

mix_data.profiles = profiles;
mix_data.proportions = proportions;
mix_data.expression=expression;
mix_data.cell_types = cell_types;
mix_data.region = region;

end


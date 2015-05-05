function [true_profiles, experiment_name, cell_type_description] =...
    load_okaty_cell_type_data(reference, region)
% Extract the expression from a specific experiment where expression was 
% extracted for each cell type seperatly.
% Since we have a limited number of astrocytes and oligo we create a mix differently:
%  -  for cortical regions we take a different sample as the true profile
%  -  for the cerebellum we take the mean profile
%  -  for the striatum we take the mean cortical astro + noise
%  -  for the brainstem we take the mean (cortical astro,cereb astro) + noise
%  -  for the spinal cord we take the mean (cortical astro,cereb astro) + noise
% 

data = load_data('okaty2011');
[num_samples, num_celltypes] = size(data.sample2type);

experiment_name = sprintf('%s_%s',reference, region);

[cell_type_ids, cell_type_sample_ids] = get_ids(experiment_name);

true_profiles = zeros(size(data.expression,1), 3);
cell_type_description = cell(length(cell_type_ids),1);
for i = 1:length(cell_type_ids)
    inds64 = ismember(data.cell_type_id, cell_type_ids{i});
    inds195 = any((data.sample2type(:, inds64)),2);
    fprintf('\t%s\n', data.cell_type_description{inds64});
    cell_type_description{i} = data.cell_type_description{inds64}; 
    
    current_expression = data.expression(:, inds195);
    if ischar(cell_type_sample_ids{i})
        switch cell_type_sample_ids{i}
            case 'mean'
                true_profiles(:, i) = mean(current_expression, 2);
            case 'mean_plus_noise'
                rand_1_minus1 = 2 * rand(size(current_expression,1),1) - 1;
                samp_var = 2 * std(current_expression,1,2) .* rand_1_minus1;
                true_profiles(:, i) = mean(current_expression, 2) + samp_var;
        end
    else
        true_profiles(:, i) = current_expression(:,cell_type_sample_ids{i});
    end
end

% Transform to linear space before adding noise
true_profiles = 2.^true_profiles;

end

function [cell_type_ids, cell_type_sample_ids] = get_ids(experiment_name)
%
    fprintf('selecting for %s:\n', experiment_name);
    switch experiment_name
      case 'Doyle_cortex_L5A'
        cell_type_ids = {...
              'H_lva_cortex',...
              'H_astro_cortex',...
              'H_mixed_olig_cortex'
              };
        cell_type_sample_ids = {1,1,1};
      case 'Doyle_cortex_L5B'
        cell_type_ids = {...
              'H_lvb_cortex',...
              'H_astro_cortex',...
              'H_mixed_olig_cortex'
              };
        cell_type_sample_ids = {2,2,2};
      case 'Doyle_cortex_L6'
        cell_type_ids = {...
              'H_l6_cortex',...
              'H_astro_cortex',...
              'H_mixed_olig_cortex'
              };
        cell_type_sample_ids = {3,3,3};
      case 'Doyle_cerebellum'
        cell_type_ids = {...
              'H_purkinje_cerebellum',...
              'H_astro_cerebellum',...
              'H_mixed_olig_cerebellum'
              };
        cell_type_sample_ids = {'mean','mean','mean'};
     case 'Doyle_striatum'
        cell_type_ids = {...
              'H_drd1+_msn_striatum',...
              'H_astro_cortex',...
              'H_mixed_olig_cortex'
              };
        cell_type_sample_ids = {1,'mean_plus_noise','mean_plus_noise'};
     case 'Doyle_spinal_cord'
        cell_type_ids = {...
              'H_chol_spinal',...
              {'H_astro_cortex','H_astro_cerebellum'},...
              {'H_mixed_olig_cortex','H_mixed_olig_cerebellum'}
              };
        cell_type_sample_ids = {1,'mean_plus_noise','mean_plus_noise'};
     case 'Doyle_brainstem'
        cell_type_ids = {...
              'H_mn_brainstem',...
              {'H_astro_cortex','H_astro_cerebellum'},...
              {'H_mixed_olig_cortex','H_mixed_olig_cerebellum'}
              };
        cell_type_sample_ids = {1,'mean_plus_noise','mean_plus_noise'};
    end

end

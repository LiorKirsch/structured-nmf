%
% Create a synthetic dataset with mix from okaty data
%

regions = {'cortex_L5A','cortex_L5B','cortex_L6',...
        'striatum','cerebellum','brainstem','spinal_cord'};

reference = 'Doyle';
mix_parms.num_tissues = 1000;
mix_parms.proportion_variability = 0.1;
mix_parms.measurement_noise = 0.01;
    
for i = 1:length(regions)
    region = regions{i};
    rng(sum(region));


    [profiles, mix_parms.experiment_name, cell_types]=...
        load_okaty_cell_type_data(reference, region);


    [proportions, expression, mix_parms.base_proportions] = mix_true_profiles(...
        profiles, region, mix_parms);


    % Save to file 
    [full_name, file_name, dir_name] = set_filenames('mixure', mix_parms);
    vars = {'profiles', 'proportions', 'expression','cell_types','region'};
    save(full_name, vars{1:end});
    fprintf('saved data to [%s/%s]\n', dir_name, file_name);
end


function rand_subset = createRandSubset(total_num_samp, num_samples_list, subsample_repeats)


rng(0); rand_subset = {};
for i_ns = 1:length(num_samples_list)
    rand_subset{i_ns} = nan(subsample_repeats, num_samples_list(i_ns) );
    for j_sr = 1:subsample_repeats
        rand_subset{i_ns}(j_sr,:) = randperm(total_num_samp,num_samples_list(i_ns));
    end
end


end
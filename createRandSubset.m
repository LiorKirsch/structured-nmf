function rand_subset = createRandSubset(total_num_samp, num_samples_list, subsample_repeats)
% Select random subsets of size 'num_samples_list' and repeat that process 
% subsample_repeat number of times.

rng(0); % to make the same selection every time.
rand_subset = cell(length(num_samples_list),1);
for i_ns = 1:length(num_samples_list)
    if length(total_num_samp) ==1
        rand_subset{i_ns} = nan(subsample_repeats, num_samples_list(i_ns) );
        for j_sr = 1:subsample_repeats
            rand_subset{i_ns}(j_sr,:) = randperm(total_num_samp,num_samples_list(i_ns));
        end
    else
        rand_subset{i_ns} = cell(length(total_num_samp),1);
        for j_total = 1:length(total_num_samp)
            rand_subset{i_ns}{j_total} = nan(subsample_repeats, num_samples_list(i_ns) );
            for j_sr = 1:subsample_repeats
                rand_subset{i_ns}{j_total}(j_sr,:) = randperm(total_num_samp(j_total),num_samples_list(i_ns));
            end 
        end
    end
end


end

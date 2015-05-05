function rand_subset = createRandSubset(total_num_samp, num_samples_list, subsample_repeats,sample_group)
% Select random subsets of size 'num_samples_list' from each group.
% It returns a cell array within each cell a matrix [time_to_repeat X num_samples*num_groups]

% if no group id is provided treat all is coming from a single group
if ~exist('sample_group','var')
    sample_group = ones(total_num_samp,1);
end
assert(length(sample_group) == total_num_samp,'each sample should have a group id');

uniq_group_id = unique(sample_group,'stable');
num_groups = length(uniq_group_id);
    
rng(0); rand_subset = cell(length(num_samples_list),1);
for i_ns = 1:length(num_samples_list)
    rand_subset{i_ns} = nan(subsample_repeats, num_samples_list(i_ns) *num_groups);
    for j_sr = 1:subsample_repeats
        curr_inds = nan(num_samples_list(i_ns), num_groups);
        
        for m_groups = 1:num_groups;
            group_inds = find(ismember(sample_group,uniq_group_id(m_groups)));
            group_num_sampl = length(group_inds);
            curr_inds(:,m_groups) = randperm(group_num_sampl,num_samples_list(i_ns));
            curr_inds(:,m_groups) = group_inds(curr_inds(:,m_groups));
        end
        rand_subset{i_ns}(j_sr,:) = curr_inds(:);
    end
end


end
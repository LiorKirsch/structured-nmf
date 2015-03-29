% Commpute the distribution of matched-corrs with random sampels
% from mix_data.expression
num_tissues = size(mix_data.expression,2);

for i=1:300
    samples = ceil(rand(1,3)*num_tissues);
    HH = mix_data.expression(:,samples)';
    [~, ~, randscores(i)] = match_profiles_to_gt(nan(size(HH)), HH, mix_data.profiles');
end
[a, x] = hist(randscores, 70);
plot(x, a);
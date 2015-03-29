%
% Create a synthetic dataset with mix from okaty data
%

data = load_data('okaty2011');

% Extract the expression from a specific experiment
experiment_name = 'Cahoy';
experiment_inds = strmatch(experiment_name, data.reference);
sample_inds = any((data.sample2type(:,experiment_inds))');
% TODO(lior) add all other fields using filter_XXX

is_neuro = double(data.sample2type) * double(data.is_neuron);
is_astro = double(data.sample2type) * double(data.is_astro);
is_oligo = double(data.sample2type) * double(data.is_oligo);

sample_neuro = find(is_neuro' .* sample_inds, 1, 'first');
sample_astro = find(is_astro' .* sample_inds, 1, 'first');
sample_oligo = find(is_oligo' .* sample_inds, 1, 'first');

samples = [sample_neuro; sample_astro; sample_oligo];
profiles = data.expression(:,samples);

% Transform to linear space before adding noise
profiles = 2.^profiles;

p_neuro = 0.6; % journal.frontiersin.org/article/10.3389/fnana.2014.00127
p_astro = 0.1;
p_oligo = 0.3;
p_all =  [p_neuro, p_astro, p_oligo];

num_tissues = 1000;
proportion_variability = 0.1;
measurement_noise = 0.01;
expression = zeros(size(data.expression,1), num_tissues);
proportions = zeros(num_tissues, 3);
rng(0);
for i_tissue = 1:num_tissues
    if mod(i_tissue,100)==0, fprintf('.');end
    p = p_all .* (randn(1,3)*proportion_variability+1);
    p = p/sum(p);
    mix = p * profiles';
    expression(:,i_tissue) = mix .* (randn(size(mix))*measurement_noise+1);
    proportions(i_tissue, 1:3) = p;
end

% Save to file 
vars = {'profiles', 'proportions', 'expression'};

dirname = fullfile('/', 'cortex', 'data', 'microarray', 'mouse', ...
                   'Okaty2011', 'Mixtures');
filename = sprintf('okaty2011-lin-lin_%s_MN%g_PR%d-%d-%d_PVAR%g.mat', ...
                   lower(experiment_name), measurement_noise, ...
                   ceil(p_all*100), proportion_variability);
save(fullfile(dirname, filename), vars{1:end});
fprintf('saved data to [%s/%s]\n', dirname, filename);


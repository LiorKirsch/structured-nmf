%
% Create a synthetic dataset with mix from okaty data
%

data = load_data('okaty2011');

% Extract the expression from a specific experiment
% experiment_name = 'Cahoy';
experiment_name = 'Doyle';
experiment_inds = strmatch(experiment_name, data.reference);
sample_inds = any((data.sample2type(:,experiment_inds))');
% TODO(lior) add all other fields using filter_XXX

is_neuro = double(data.sample2type) * double(data.is_neuron);
is_astro = double(data.sample2type) * double(data.is_astro);
is_oligo = double(data.sample2type) * double(data.is_oligo);

sample_selection_for_mix = take_from_struct(parms, 'sample_selection_for_mix', '');
if isempty(sample_selection_for_mix)
    allinds = ones(size(is_neuro));
    suffix = '';
else
    switch experiment_name
      case 'Doyle'
        inds64 = strmatch('H_cck+_cortex', data.cell_type_id);
        allinds = any((data.sample2type(:, inds64)),2);
        
        inds64 = strmatch('H_astro_cortex', data.cell_type_id);
        allinds = allinds | any((data.sample2type(:, inds64)),2);

        inds64 = strmatch('H_mixed_olig_cortex', data.cell_type_id);
        allinds = allinds | any((data.sample2type(:, inds64)),2);
        suffix = '-mixedcortex' 
    end
end

% Create teh profiles
samples_neuro = find(is_neuro' .* sample_inds .* allinds');
samples_astro = find(is_astro' .* sample_inds .* allinds');
samples_oligo = find(is_oligo' .* sample_inds .* allinds');

profiles  = zeros(size(data.expression,1), 3);
profiles(:, 1) = mean(data.expression(:, samples_neuro), 2);
profiles(:, 2) = mean(data.expression(:, samples_astro), 2);    
profiles(:, 3) = mean(data.expression(:, samples_oligo), 2);


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
filename = sprintf('okaty2011-lin-lin_%s%s_MN%g_PR%d-%d-%d_PVAR%g.mat', ...
                   lower(experiment_name), lower(suffix), measurement_noise, ...
                   ceil(p_all*100), proportion_variability);
save(fullfile(dirname, filename), vars{1:end});
fprintf('saved data to [%s/%s]\n', dirname, filename);


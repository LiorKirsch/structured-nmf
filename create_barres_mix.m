%
% Create a synthetic dataset with mix from barres data
%

data = load_data('barres2014');

% the data is already in linear scale

% cell_types = {'Astrocytes','Neuron','Oligodendrocyte Precursor Cell','Newly Formed Oligodendrocyte','Myelinating Oligodendrocytes','Microglia','Endothelial Cells'};
profiles = [data.data(:,strcmp('Neuron',data.cell_types)) ,...
            data.data(:,strcmp('Astrocytes',data.cell_types)) ,...
            data.data(:,strcmp('Myelinating Oligodendrocytes',data.cell_types))];

p_neuro = 0.6; % journal.frontiersin.org/article/10.3389/fnana.2014.00127
p_astro = 0.1;
p_oligo = 0.3;
p_all =  [p_neuro, p_astro, p_oligo];

num_tissues = 1000;
proportion_variability = 0.1;
measurement_noise = 0.1;
expression = zeros(size(data.data,1), num_tissues);
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

dirname = fullfile('/', 'cortex', 'data', 'RNA-Seq', 'mouse', ...
                         'Barres-2014', 'Mixtures');
filename = sprintf('barres2014-lin-lin_MN%g_PR%d-%d-%d_PVAR%g.mat', ...
                   measurement_noise, ...
                   ceil(p_all*100), proportion_variability);
save(fullfile(dirname, filename), vars{1:end});
fprintf('saved data to [%s/%s]\n', dirname, filename);


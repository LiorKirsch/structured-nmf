function priors = load_priors(dataset_name, cell_types, parms)
%
% Take priors from Okaty.
%

  
  data = load_data('okaty2011');
  if any(isnan(data.expression(:)))
      error('Data contains NaN');
  end
  
  % Extract the expression from a specific experiment
  experiment_inds64 = strmatch(dataset_name, data.reference);
  if isempty(experiment_inds64)
      fprintf('invalid dataset_name [%s]\n', dataset_name);
      fprintf('valid datasets are:\n');
      disp(unique(data.reference));
      error('bye bye');
  end
  experiment_inds195 = any((data.sample2type(:,experiment_inds64))');

  is_neuro = double(data.sample2type) * double(data.is_neuron);
  is_astro = double(data.sample2type) * double(data.is_astro);
  is_oligo = double(data.sample2type) * double(data.is_oligo);

  prior_selection = take_from_struct(parms, 'prior_selection', ...
                                     'mean_cortex');  
  fprintf('load_priors: use %s to create priors\n', prior_selection);
  switch prior_selection
    case 'mean',
        sample_neuro = logical(is_neuro' .* sample_inds);
        sample_astro = logical(is_astro' .* sample_inds);
        sample_oligo = logical(is_oligo' .* sample_inds);
    
    case 'mean_cortex',
        region_inds64 = strmatch('cortex', lower(data.anatomical_region));
        region_inds195 = any((data.sample2type(:,region_inds64))');

        sample_neuro = logical(is_neuro' .* experiment_inds195 .* region_inds195);
        sample_astro = logical(is_astro' .* experiment_inds195 .* region_inds195);
        sample_oligo = logical(is_oligo' .* experiment_inds195 .* region_inds195);    
    
    case 'mean_cerebellum',
        region_inds64 = strmatch('cerebellum', lower(data.anatomical_region));
        region_inds195 = any((data.sample2type(:,region_inds64))');

        sample_neuro = logical(is_neuro' .* experiment_inds195 .* region_inds195);
        sample_astro = logical(is_astro' .* experiment_inds195 .* region_inds195);
        sample_oligo = logical(is_oligo' .* experiment_inds195 .* region_inds195);
        
    case 'single'
      error('Not supported yet');
  end
  
  priors = zeros(size(data.expression, 1),0);
  
  for i=1:length(cell_types)
      switch cell_types{i}
        case {'neuron', 'neuro'}, 
          priors(:, end+1) = mean(data.expression(:,sample_neuro), 2);
          fprintf('averaged %d samples for neurons\n', nnz(sample_neuro));
          assert(nnz(sample_neuro)>0, 'No neuron samples found');        
        case 'astro',
          priors(:, end+1) = mean(data.expression(:,sample_astro),2);          
          fprintf('averaged %d samples for astro\n', nnz(sample_astro));
          assert(nnz(sample_astro)>0, 'No astro samples found');
        case 'oligo', 
          priors(:, end+1) = mean(data.expression(:,sample_oligo),2);          
          fprintf('averaged %d samples for oligo\n', nnz(sample_oligo));
          assert(nnz(sample_oligo)>0, 'No oligo samples found');
        otherwise, 
          error('invalid cell type = [%s]\n', cell_types{i});
      end          
  end
end

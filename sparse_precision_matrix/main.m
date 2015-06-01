addpath('/cortex/code/cellmix/');  

init;

% parms = conf(parms);

dataset = take_from_struct(parms, 'dataset', 'brainspan2014');

dataset = 'zapala2005';
[default_regions, parms.species] = get_region_set(dataset);
regions = take_from_struct(parms, 'regions', default_regions);
regions = sort(regions);

switch dataset
    case 'kang2011',      %==== Kang ===
      parms.dataset_file = 'kang_regions';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('kangCortexAndStriatum', []);
      
    case 'brainspan2014', %==== Brainspan ===
      parms.dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions,'_'));
      [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
          = load_expression_and_regions('brainspan_rnaseq', regions);
    
    case 'zapala2005',    %==== Zapala selected regions ===
      parms.dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
      [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
       ~] = load_expression_and_regions('zapalaMouse', regions);
      gross_structures_info{strcmp(gross_structures_info, 'Isocortex')} ...
          = 'Cerebral_cortex';

      case 'human6' , %==== Human6 selected regions ===
        parms.dataset_file = 'Human6_selected_regions';
        [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
         ~] = load_expression_and_regions('human6LimitRegions', regions);
        gross_structures_info = gross_structures_info(:,4);
        gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x), ...
                                        gene_info.entrez_ids, ...
                                        'UniformOutput',false);
    otherwise 
      error('invalid dataset = [%s]\n', dataset);
end


A = cov(expression);
precision_lambda = 1;
rho = 1;
positive = 1;
[precision_matrix, Y, iter, positivity_fail]=inverse_cov(A,precision_lambda,rho,positive);

[Theta, W] = graphicalLasso(A, precision_lambda, 1000, 10^-7);

file_name = set_filenames('precision',parms);
save(file_name, 'precision_matrix','gene_info');
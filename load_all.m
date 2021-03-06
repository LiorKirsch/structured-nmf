function [dataset_file, expression, gross_region_vec, gene_info, ...
          gross_structures_info] = load_all(dataset, regions,just_file_name)

      if ~exist('just_file_name','var')
          just_file_name = false;
      end
      
      if just_file_name
        
        switch dataset
            case 'kang2011',      %==== Kang ===
              dataset_file = 'kang_regions';
            case 'brainspan2014', %==== Brainspan ===
              dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions,'_'));
            case 'zapala2005',    %==== Zapala selected regions ===
              dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
            case 'human6' , %==== Human6 selected regions ===
              dataset_file = 'Human6_selected_regions';
            otherwise 
              error('invalid dataset = [%s]\n', dataset);
        end
        expression = nan;
        gross_region_vec = nan;
        gene_info = nan;
        gross_structures_info = nan;
      else
          
        switch dataset
            case 'kang2011',      %==== Kang ===
              dataset_file = 'kang_regions';
              [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
               ~] = load_expression_and_regions('kangCortexAndStriatum', []);

            case 'brainspan2014', %==== Brainspan ===
              dataset_file = sprintf('brainspan_rnaseq_%s', strjoin(regions,'_'));
              [expression, gross_region_vec, gene_info, ~, gross_structures_info] ...
                  = load_expression_and_regions('brainspan_rnaseq', regions);

            case 'zapala2005',    %==== Zapala selected regions ===
              dataset_file = 'Zapala_isocortex_medulla_striatum_cerebellum';
              [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
               ~] = load_expression_and_regions('zapalaMouse', regions);
              gross_structures_info{strcmp(gross_structures_info, 'Isocortex')} ...
                  = 'Cerebral_cortex';

            case 'human6' , %==== Human6 selected regions ===
              dataset_file = 'Human6_selected_regions';
              [expression, gross_region_vec, gene_info, ~, gross_structures_info, ...
               ~] = load_expression_and_regions('human6LimitRegions', regions);
              gross_structures_info = gross_structures_info(:,4);
              gene_info.entrez_ids = arrayfun(@(x) sprintf('%d',x), ...
                                              gene_info.entrez_ids, ...
                                              'UniformOutput',false);
            otherwise 
              error('invalid dataset = [%s]\n', dataset);
        end

      end
end
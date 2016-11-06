function data = load_data(dataset_name, parms)
%
%
%    
    persistent local_data
    persistent local_datasetname
    if ~isempty(local_data) && strcmp(local_datasetname, ...
                                      dataset_name)
        data = local_data;
        return
    end
    
    
  switch dataset_name
    case 'okaty2011',
      dirname = fullfile('/', 'cortex', 'data', 'microarray', 'mouse', ...
                         'Okaty2011');
      filename = 'mouse_cell_type_profiles.mat';
      data = load(fullfile(dirname, filename));
      
    case 'barres2014', 
      dirname = fullfile('/', 'cortex', 'data', 'RNA-Seq', 'mouse', ...
                         'Barres-2014');
      filename = 'barres_rnaseq.mat';
      data = load(fullfile(dirname, filename));
      
    case {'okaty2011_cahoy_MN0.1_PR60-10-30_PVAR0.1',...
          'okaty2011-lin-lin_cahoy_MN0.01_PR60-10-30_PVAR0.1',...
          'okaty2011-lin-lin_cahoy_MN0.05_PR60-10-30_PVAR0.1',...
          'okaty2011-lin-lin_cahoy_MN0.1_PR60-10-30_PVAR0.1'}
      dirname = fullfile('/', 'cortex', 'data', 'microarray', 'mouse', ...
                         'Okaty2011', 'Mixtures');
      mixname = strrep(dataset_name, 'okaty2011_', '');
      filename = sprintf('%s.mat', dataset_name);
      data = load(fullfile(dirname, filename));      
      
    case {'barres2014-lin-lin_MN0.01_PR60-10-30_PVAR0.1',...
          'barres2014-lin-lin_MN0.1_PR60-10-30_PVAR0.1' }
      dirname = '/cortex/data/RNA-Seq/mouse/Barres-2014/Mixtures/';
%       mixname = strrep(dataset_name, 'barres2014-', '');
      filename = sprintf('%s.mat', dataset_name);
      data = load(fullfile(dirname, filename));    
      
    case {'darmanis2015'}
      dirname = '/cortex/data/RNA-Seq/human/Darmanis2015/';
      filename = 'rnaseq_celltypes_GPL18573.mat';
      data = load(fullfile(dirname, filename));
      
      data.expression_matrix = data.expression_matrix';
      data.expression = data.expression_matrix;
      data.all_symbols = data.gene_names;
      data.gene_symbol = data.gene_names;
      data.refer_to_index = 1:length(data.gene_symbol);
              
      
    case 'kang2011', 
      error('not supported yet\n');      
    
    
      
    otherwise
      error('invalid dataset_name [%s]\n', dataset_name);
  end
  
  local_data = data;
  local_datasetname = dataset_name;
  
end
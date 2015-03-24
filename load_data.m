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
      
    case 'okaty2011_cahoy_MN0.1_PR60-10-30_PVAR0.1',
      dirname = fullfile('/', 'cortex', 'data', 'microarray', 'mouse', ...
                         'Okaty2011', 'Mixtures');
      mixname = strrep(dataset_name, 'okaty2011_', '');
      filename = sprintf('%s.mat', dataset_name);
      data = load(fullfile(dirname, filename));      
      
      
    case 'kang2011', 
      error('not supported yet\n');      
    
    
    otherwise
      error('invalid dataset_name [%s]\n', dataset_name);
  end
  
  local_data = data;
  local_datasetname = dataset_name;
  
end
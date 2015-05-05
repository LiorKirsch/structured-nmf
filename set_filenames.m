function [full_name, file_name, dir_name] = set_filenames(file_type, parms)
%
  switch file_type
      
    case 'demixing', 
      dir_name = fullfile('/cortex/users/lior/nmf/runs/',parms.dataset_file );
      [parmstr, dirparmstr] = set_parmstr(parms);
      dir_name = fullfile(dir_name,dirparmstr);
      file_name = sprintf('demix_%s.mat', parmstr);
      auto_mkdir = 1;
    case 'figure'
        dir_name = fullfile('figures',parms.dataset_file);
        parmstr = set_parmstr(parms);
        file_name = sprintf('%s_%s_%s.png', parms.fig_x_axis, parms.corr_type,parmstr);
        auto_mkdir = 1;
    case 'proportions-figure'
        dir_name = fullfile('figures',parms.dataset_file);
        parmstr = set_parmstr(parms);
        file_name = sprintf('proportions_%s_%s_%s.png', parms.fig_x_axis, parms.corr_type ,parmstr);
        auto_mkdir = 1;
    case 'mixure'    
        dir_name = fullfile('/', 'cortex', 'data', 'microarray', 'mouse', ...
                   'Okaty2011', 'Mixtures');
               
        file_name = sprintf('okaty2011-%s_MN%g_PR%d-%d-%d_PVAR%g.mat', ...
                   lower(parms.experiment_name), parms.measurement_noise, ...
                   ceil(parms.base_proportions*100), parms.proportion_variability);
        auto_mkdir = 1;
    otherwise 
      error('invalid file_type = [%s]\n', file_type);
  end


  if exist(dir_name, 'dir') ~= 7
    if auto_mkdir
      fprintf('\nset_filenames.m:\n');    
      fprintf('\tAuto-Make directory "%s" for filename %s\n',dir_name, ...
	      file_name);
      mkdir(dir_name);
    end
  end


  full_name  = fullfile(dir_name, file_name);
end




function [full_name, file_name, dir_name] = set_filenames(file_type, parms)
%
  switch file_type
      
    case 'demixing', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/runs/');
      [parmstr, dirparmstr] = set_parmstr(parms);
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('demix_%s.mat', parmstr);
      auto_mkdir = 1;
    case 'baselines', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/baselines/');
      [parmstr, dirparmstr] = set_parmstr(parms);
      if  isfield(parms, 'true_dataset')
          if ~strcmp(parms.true_dataset, 'okaty')
              parmstr = sprintf('%s_%s', parms.true_dataset, parmstr);
          end
      end
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('baseline_%s.mat', parmstr);
      auto_mkdir = 1;
    
    case 'demixing_rand_restart', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/runs/rand_restarts/');
      [parmstr, dirparmstr] = set_parmstr(parms);
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('restart_%s.mat', parmstr);
      auto_mkdir = 1;
    case 'results', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/results/');
      [parmstr, dirparmstr] = set_parmstr(parms);
      if  isfield(parms, 'true_dataset')
          if ~strcmp(parms.true_dataset, 'okaty')
              parmstr = sprintf('%s_%s', parms.true_dataset, parmstr);
          end
      end
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('results_%s.mat', parmstr);
      auto_mkdir = 1;
    case 'precision', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/precision/');
      [parmstr, dirparmstr] = set_parmstr(parms);
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('precision_%s.mat', parmstr);
      auto_mkdir = 1;
    case 'gene_subset', 
      base_dir_name = fullfile('/cortex/users/lior/nmf/runs/gene_subset');
      [parmstr, dirparmstr] = set_parmstr(parms);
      dir_name = fullfile(base_dir_name, dirparmstr);
      file_name = sprintf('gene_subset%s.mat', parmstr);
      auto_mkdir = 1;
    case 'figure'
        base_dir_name = 'figures';
        [parmstr, dirparmstr] = set_parmstr(parms);
        dir_name = fullfile(base_dir_name, dirparmstr);
        file_name = sprintf('%s_%s_%s.png', parms.fig_x_axis, parms.corr_type,parmstr);
        auto_mkdir = 1;
    case 'figure_real'
        base_dir_name = 'figures_real';
        [parmstr, dirparmstr] = set_parmstr(parms);
        dir_name = fullfile(base_dir_name, dirparmstr);
        file_name = sprintf('%s_%s_%s.png', parms.fig_x_axis, parms.corr_type,parmstr);
        auto_mkdir = 1;
    case 'proportions-figure'
        dir_name = fullfile('figures',parms.dataset_file);
        parmstr = set_parmstr(parms);
        file_name = sprintf('proportions_%s_%s_%s.png', parms.fig_x_axis, parms.corr_type ,parmstr);
        auto_mkdir = 1;
     case 'figure_confusion'
        [parmstr, dirparmstr] = set_parmstr(parms);
%         dir_name = fullfile('figures',parms.dataset_file);
        file_name = sprintf('confusion_%s.png',parmstr);
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
      fprintf('\tAuto-Make directory "%s" for filename %s\n', dir_name, ...
              file_name);
      mkdir(dir_name);
    end
  end


  full_name  = fullfile(dir_name, file_name);
end




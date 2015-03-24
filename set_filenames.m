function [full_name, file_name, dir_name] = set_filenames(file_type, parms)
%
  switch file_type
      
    case 'demixing', 
      dir_name = '/cortex/users/lior/nmf/runs/';
      parmstr = set_parmstr(parms);
      file_name = sprintf('demix_%s.mat', parmstr);
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




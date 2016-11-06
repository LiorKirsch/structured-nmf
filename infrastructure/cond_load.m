
function [do_calc, varargout] = cond_load(filename, do_force, varargin)
%
% Load results from file if possible & allowed. 
%
% [do_calc, stats] = conditional_load(filename, do_force, varargin)
%
% Example: 
%   vars  = {'x','y'}
%   [do_calc, x, y] = cond_load(filename, do_force, vars{1:end});
%   if(do_calc < 1 ), return; end
%
% Example2: 
%   vars  = {'x','y', 'quiet'}
%   [do_calc, x, y] = cond_load(filename, do_force, vars{1:end});
%
%


% Debug
%  disp(filename)
%  dbstack
  
  % Special verbosity signal
  if isnan(do_force)
      verbosity = false;
  else
      verbosity = true;
  end

  is_quiet = cellfun(@(x) strcmp(x,'quiet'), varargin);
  if any(is_quiet)
      varargin = varargin(not(is_quiet));
      verbosity = false;
  end

  if do_force>0
    do_calc=1;
    n_vars = length(varargin);    
    for i_var = 1: n_vars
      var = varargin{i_var};
      cmd = sprintf('   varargout{%d}=[];',i_var,var);    
      eval(cmd);
    end
    return
  end  
  
  % Init var
  vars = varargin;
  n_vars = length(vars);
  for i_var = 1: n_vars 
    cmd = sprintf('%s = [];',vars{i_var}); eval(cmd);
  end
  varargout=[];

  if nargout-1 ~= length(vars)
    fprintf('Warning: arguments mismatch in cond_load.m\n');
    fprintf('Calling stack is:\n');
    dbstack
    fprintf('End of stack\n');    
    error('arguments mismatch in cond_load.m');
  end
  
  % Check if file exists
  do_calc = ~logical(exist(filename,'file'));

  S.dummy=[];
  if do_calc<1
      try
	S = load(filename);
      catch
	fprintf('filename = "%s"\n',filename);
	fprintf('Warning: File exist, but could not be read\n');
	do_calc=1;
	keyboard	
      end
  end

  if exist(filename, 'file') == 0
    if verbosity
      filename_prt = shorten_filename(filename);
      fprintf('File "%s" doesnt exist.\n',filename_prt);
      fprintf('> Calculate all vars %s\n',datestr(now,13));
    end
    for i_var = 1: n_vars
      var = vars{i_var};
      cmd = sprintf('   varargout{%d}=[];',i_var,var);    
      eval(cmd);
    end
    return
  end

  % File exists, Check if var was read
  for i_var = 1: n_vars
    var = vars{i_var};
    if isfield(S, var)==0
      do_calc=1;      
      if verbosity
          fprintf('No variable [%s] found. Recalculate\n', var);
      end
      varargout{i_var} = [];
    else
      varargout{i_var} = S.(var);
    end
  end  
end
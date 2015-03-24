function short = shorten_filename(long)
%
% short = shorten_filename(long)
%
% Replace strings in paths for printing short file names
% (C) Gal Chechik 2007 
% (C) Gal Chechik 2013 
% 
  short = long; 
  % short = strrep(short,'/home/gal/','~/');
  short = strrep(short,'/home/lab/gal/Projects/','//');  
  short = strrep(short,'/home/lab/gal/','~/');  
  % short = strrep(short,'/afs/cs/u/gal/','~/');
  % short = strrep(short,'/afs/cs.stanford.edu/u/gal/','~/');
  % short = strrep(short,'/Projects/TimeCourses/Runs/Sources_M','/TC');
  % short = strrep(short,'/Projects/TimeCourses/Runs/SingleGenesFit/Imodel/Pathways/R14_Num_NOCPLX/','/TCR/');
  % short = strrep(short,'/Projects/TimeCourses/Runs','/TC/..');
  % short = strrep(short,'/Data/MetPaths/Yeast/Palsson2003/','/Palsson/');  

end

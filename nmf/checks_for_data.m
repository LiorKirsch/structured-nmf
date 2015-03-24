function checks_for_data(data_matrix)


assert(all(data_matrix(:) >=0),'Error, Data matrix values should be geq than zero');
assert(~any(sum(data_matrix,1)==0),'Data matrix should not contain a zero column');
assert(~any(sum(data_matrix,2)==0),'Data matrix should not contain a zero row');

  

end
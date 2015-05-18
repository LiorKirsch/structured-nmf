function new_expression_target = match_expression_across_experiments(expression_base, expression_target)
% Changes the values of "expression_target" in a way that they will match
% the distribution of values in expression_base
% It replaces the lowest values of target with the lowest value of base
% Then, it replaces the second lowest value of target ...
% Until all the values in target are the values in base.
%
% example:
%     match_expression_across_experiments([2 4 6 8 10 12], [14 12 18; 20 1 3])
%      -->>
%     ans =     8     6    10
%              12     2     4


    [sorted_exp_A, inds_A] = sort( expression_base(:) );
    [sorted_exp_B, inds_B] = sort( expression_target(:) );
    
    [~,reverse_inds_B] = sort(inds_B);
    new_expression_target = sorted_exp_A(reverse_inds_B);
    
    new_expression_target = reshape(new_expression_target, size(expression_target));

end
function [reorder_A, reorder_B] = reorderUsingId(IdA, IdB)
    intersection_ids = intersect(IdA, IdB);
    [~, reorder_A] = ismember(intersection_ids, IdA);
    [~, reorder_B] = ismember(intersection_ids, IdB);
end
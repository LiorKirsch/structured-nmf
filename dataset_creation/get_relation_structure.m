function [relationMatrix,all_strct] = get_relation_structure(limit_to,relation_type)

    if ~exist('relation_type','var')
        relation_type = 'common_parent_level';
    end
    
    [tree_structure_matrix,all_strct] = get_tree_structure(false,limit_to);
    tree_structure_matrix = sparse(double(tree_structure_matrix));

    switch relation_type
        case 'relations_dist'
            [~,relationMatrix] = computeDistanceBetweenNodes(tree_structure_matrix);
            relationMatrix = exp(-1*relationMatrix);
        case 'relations_parentdist'
            directedDistanceMatrix = computeDistanceBetweenNodes(tree_structure_matrix);
            [ closestCommonParentIndex, meanDistanceToParent, shortDistanceToParent, longDistanceToParent] = ...
                distanceToCommonParent(tree_structure_matrix, directedDistanceMatrix);
            relationMatrix = longDistanceToParent;
            relationMatrix = exp(-1*relationMatrix);
        case 'relations_parent_level'
            [~, nodeLevel] = allChildNodes(tree_structure_matrix);
            directedDistanceMatrix = computeDistanceBetweenNodes(tree_structure_matrix);
            closestCommonParentIndex = distanceToCommonParent(...
                tree_structure_matrix, directedDistanceMatrix);
            
            closestCommonParentIndex(1:size(closestCommonParentIndex,1)+1:end) = 1;
            relationMatrix = nodeLevel(closestCommonParentIndex);
            relationMatrix(1:size(relationMatrix,1)+1:end) = 0;
            
            relationMatrix = exp(relationMatrix);
        otherwise
            error('unkown relation type');
    end


    mask = ismember(all_strct, limit_to);
    
    relationMatrix = relationMatrix(mask, mask);
    all_strct = all_strct(mask);

%     figure;imagesc(relationMatrix);colorbar; colormap(jet);
%     ax = gca;
%     ax.XTickLabel = all_strct;
%     ax.XTickLabelRotation	=45;
%     
    

end
function [directedDistanceMatrix,unDirectedDistanceMatrix] =...
        computeDistanceBetweenNodes(dependecyMatrix)
% This function computes distances between nodes in the graph
% It uses the matlab_bgl library.

    addpath('~/Projects/matlab_bgl')
    undirectedMatrix = dependecyMatrix + dependecyMatrix';
    directedDistanceMatrix = nan(size(dependecyMatrix));
    unDirectedDistanceMatrix = nan(size(dependecyMatrix));
    for i = 1:size(dependecyMatrix,1)
        [nodeDistance, ~] = dijkstra_sp(dependecyMatrix,i);
        directedDistanceMatrix(:,i) = nodeDistance;
        
        [nodeDistance ~] = dijkstra_sp(undirectedMatrix,i);
        unDirectedDistanceMatrix(:,i) = nodeDistance;
        
    end
end

function [allChilds, nodeLevel] = allChildNodes(dependencyMatrix)
    allChilds = inv(eye(size(dependencyMatrix)) - dependencyMatrix);
    nodeLevel = sum(allChilds,1);
end

function [ closestCommonParentIndex, meanDistanceToParent, shortDistanceToParent, longDistanceToParent] = ...
    distanceToCommonParent(dependencyMatrix, directedDistanceMatrix)
    [allChilds, nodeLevel] = allChildNodes(dependencyMatrix);

    

    closestCommonParentIndex = zeros(size(directedDistanceMatrix));
    meanDistanceToParent = zeros(size(directedDistanceMatrix));
    shortDistanceToParent = zeros(size(directedDistanceMatrix));
    longDistanceToParent = zeros(size(directedDistanceMatrix));

    numberOfNodes = size(directedDistanceMatrix,1);


    for i = 1: numberOfNodes
        for j= i+1 : numberOfNodes
            ancestorOf_i_j = [ allChilds(:,i) , allChilds(:,j)];
            ancestorOf_i_j = all( ancestorOf_i_j,2);
            ancestorsIndices = find(ancestorOf_i_j);

            [~, bestParentmInIndex] = max(nodeLevel (ancestorOf_i_j) );
            bestParentmInIndex = ancestorsIndices(bestParentmInIndex);

            parentDistanceFrom_i = directedDistanceMatrix(i,bestParentmInIndex);
            parentDistanceFrom_j = directedDistanceMatrix(j,bestParentmInIndex);
            distances = [parentDistanceFrom_i, parentDistanceFrom_j];
            sortedDist = sort(distances);
            shortDistanceToParent(i,j) = sortedDist(1);
            longDistanceToParent(i,j) = sortedDist(2);
            meanDistanceToParent(i,j) = mean(distances);
            closestCommonParentIndex(i,j) = bestParentmInIndex;

        end
    end

    closestCommonParentIndex = closestCommonParentIndex + closestCommonParentIndex';
    meanDistanceToParent = meanDistanceToParent + meanDistanceToParent';
    shortDistanceToParent = shortDistanceToParent + shortDistanceToParent';
    longDistanceToParent = longDistanceToParent + longDistanceToParent';

end
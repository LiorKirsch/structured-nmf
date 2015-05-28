function [relationMatrix, all_strct] = ...
        get_relation_structure(tree_structure_matrix, all_strct, ...
                                             limit_to, relation_type, X)
% get the relation matrix
% the larger the weight the strong the connection between the region is
% weight of zero means no connection between the regions
%
    if ~exist('relation_type','var')
        relation_type = 'common_parent_level';
    end
    
%   [tree_structure_matrix, all_strct] = get_tree_structure(false,limit_to);
    tree_structure_matrix = sparse(double(tree_structure_matrix));

    switch relation_type
        case 'tree'
            only_child = find(sum(tree_structure_matrix, 2) ==1);
            for i =1:length(only_child)
                curr_ind = only_child(i);
                child_ind = find(tree_structure_matrix( curr_ind, :));
                parent_ind = find(tree_structure_matrix( :, curr_ind));
                tree_structure_matrix(parent_ind, child_ind) = true;
            end
            tree_structure_matrix(only_child,:) = [];
            tree_structure_matrix(:,only_child) = [];
            all_strct(only_child) = [];
            relationMatrix = full(tree_structure_matrix);

            disp('using tree structure');

        case 'relations_dist'
            [~, relationMatrix] = computeDistanceBetweenNodes(tree_structure_matrix);
            relationMatrix = exp(-1*relationMatrix);
            relationMatrix(1:size(relationMatrix,1)+1:end) = 0;  % the diagonal should be zero

        case 'relations_parentdist'
            directedDistanceMatrix = computeDistanceBetweenNodes(tree_structure_matrix);
            [ closestCommonParentIndex, meanDistanceToParent, ...
              shortDistanceToParent, longDistanceToParent] = ...
                distanceToCommonParent(tree_structure_matrix, ...
                                       directedDistanceMatrix);
            relationMatrix = longDistanceToParent;
            relationMatrix = exp(-1*relationMatrix);
            relationMatrix(1:size(relationMatrix,1)+1:end) = 0;  % the diagonal should be zero

        case 'relations_parent_level'
            [~, nodeLevel] = allChildNodes(tree_structure_matrix);
            directedDistanceMatrix = computeDistanceBetweenNodes(tree_structure_matrix);
            closestCommonParentIndex = distanceToCommonParent(...
                tree_structure_matrix, directedDistanceMatrix);
            
            closestCommonParentIndex(1:size(closestCommonParentIndex,1)+1:end) = 1;
            relationMatrix = nodeLevel(closestCommonParentIndex);
            
            relationMatrix = exp(relationMatrix);
            relationMatrix(1:size(relationMatrix,1)+1:end) = 0;  % the diagonal should be zero

         case 'relations_on_expression'
            relationMatrix = nan(length(X),length(X));
            for i =1:length(X)
                for j=1:length(X)
                relationMatrix(i,j) = corr(mean(X{i},1)' ,mean(X{j},1)' );
                relationMatrix(j,i) = corr(mean(X{i},1)' ,mean(X{j},1)' );
                end
            end
            relationMatrix(1:size(relationMatrix,1)+1:end) = 0;  % the diagonal should be zero
        otherwise
            error('unkown relation type');
    end


    if ~strcmp('tree', relation_type)
        mask = ismember(all_strct, limit_to);

        relationMatrix = relationMatrix(mask, mask);
        all_strct = all_strct(mask);
    end

%     figure;imagesc(relationMatrix);colorbar; colormap(jet);
%     ax = gca;
%     ax.XTick = 1:length(all_strct);
%     ax.XTickLabel = all_strct;
%     ax.XTickLabelRotation	=45;
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
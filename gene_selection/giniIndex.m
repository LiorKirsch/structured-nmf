function [ bestSplitPoints,giniAttrs,sortedGinis,sortedGiniAttrs ] = giniIndex( X,Y )
% USING THE ATTRIBUTE SELECTION MEASURE GINI INDEX
% Calculating the expected information needed to classify a tuple in data.
% https://github.com/iraquitan/matlab-scripts/tree/master/Decision%20Tree
classes = unique(Y);
nsamples = length(Y);
giniData = 0;
for i=1:length(classes)
    prob_i = sum(ismember(Y, classes{i})) / nsamples;
    giniData = giniData - prob_i^2;
end
giniData = 1 + giniData;
% Calculating how much more information still needed after partitioning the
% data to arrive an exact classification.
bestSplitPoints = cell(nsamples-1,1);
giniAttrs = zeros(nsamples-1,1);
[~,numAttr] = size(X);
for attr=1:numAttr
    [A,index] = sortrows(X,attr);
    Y_sorted = Y(index);
    midPoints = cell(nsamples-1,1);
    giniAttrDatas = zeros(nsamples-1,1);    
    for i=1:(nsamples-1)
        giniAttrData = 0;
        midPoint = (A(i,attr) + A(i+1,attr)) / 2;
        midPoints{i,1} = midPoint;
        for j=1:2
            if j == 1
                dataj = find(A(:,attr)<=midPoint);
            elseif j == 2
                dataj = find(A(:,attr)>midPoint);
            end
            mod_dataj = length(dataj);
            giniDataj = 0;
            classesj = unique(Y_sorted(dataj));
            for m=1:length(classesj)
                prob_m = sum(ismember(Y_sorted(dataj), classesj{m})) / mod_dataj;
                giniDataj = giniDataj - prob_m^2;
            end
            giniDataj = 1 + giniDataj;
            giniAttrData = giniAttrData + ((mod_dataj / nsamples) * giniDataj);
        end
        giniAttrDatas(i) = giniAttrData;
    end
    minGiniAttrData = min(giniAttrDatas);
    bestSplitPoints{attr} = midPoints{giniAttrDatas == min(giniAttrDatas)};
    % Calculating how much would be gained by branching on attribute
    giniAttrs(attr) = giniData - minGiniAttrData;
end
[sortedGinis,sortedGiniAttrs] = sortrows(giniAttrs,-1);

end
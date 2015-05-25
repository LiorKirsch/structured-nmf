function [ bestSplitPoints,gainRatioAttrs,sortedGainRatios,sortedGainRatioAttrs ] = gainRatio( X,Y,gain )
% USING THE ATTRIBUTE SELECTION MEASURE GAIN RATIO
% Calculating the expected information needed to classify a tuple in data.
% https://github.com/iraquitan/matlab-scripts/tree/master/Decision%20Tree
nsamples = length(Y);
bestSplitPoints = cell(nsamples-1,1);
gainRatioAttrs = zeros(nsamples-1,1);
[~,numAttr] = size(X);
for attr=1:numAttr
    printPercentCounter(attr, numAttr);
    A = sortrows(X,attr);
    midPoints = cell(nsamples-1,1);
    splitInfoAttrDatas = zeros(nsamples-1,1);    
    for i=1:(nsamples-1)
        splitInfoAttrData = 0;
        midPoint = (A(i,attr) + A(i+1,attr)) / 2;
        midPoints{i,1} = midPoint;
        for j=1:2
            if j == 1
                dataj = find(A(:,attr)<=midPoint);
            elseif j == 2
                dataj = find(A(:,attr)>midPoint);
            end
            mod_dataj = length(dataj);
            splitInfoAttrData = splitInfoAttrData - ((mod_dataj / nsamples) * log2(mod_dataj / nsamples));
        end
        splitInfoAttrDatas(i) = splitInfoAttrData;
    end
    minSplitInfoAttrData = min(splitInfoAttrDatas);
    bestSplitPoints{attr} = midPoints{splitInfoAttrDatas == min(splitInfoAttrDatas)};
    % Calculating how much would be gained by branching on attribute
    gainRatioAttrs(attr) = gain(attr) / minSplitInfoAttrData;
end
[sortedGainRatios,sortedGainRatioAttrs] = sortrows(gainRatioAttrs,-1);
end
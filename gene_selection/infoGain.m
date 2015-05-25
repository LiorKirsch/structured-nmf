function [bestSplitPoints, gainAttrs, sortedGains, sortedAttr] = infoGain( X,Y )
% USING THE ATTRIBUTE SELECTION MEASURE INFORMATION GAIN
% Calculating the expected information needed to classify a tuple in data.
% https://github.com/iraquitan/matlab-scripts/tree/master/Decision%20Tree
classes = unique(Y);
nsamples = length(Y);
infoData = 0;
for i=1:length(classes)
    prob_i = sum(ismember(Y, classes{i})) / nsamples;
    infoData = infoData - prob_i * log2(prob_i);
end

% Calculating how much more information still needed after partitioning the
% data to arrive an exact classification.
bestSplitPoints = cell(nsamples-1,1);
gainAttrs = zeros(nsamples-1,1);
[~,numAttr] = size(X);
for attr=1:numAttr
    printPercentCounter(attr, numAttr);
    [A,index] = sortrows(X,attr);
    Y_sorted = Y(index);
    midPoints = cell(nsamples-1,1);
    infoAttrDatas = zeros(nsamples-1,1);    
    for i=1:(nsamples-1)
        infoAttrData = 0;
        midPoint = (A(i,attr) + A(i+1,attr)) / 2;
        midPoints{i,1} = midPoint;
        for j=1:2
            if j == 1
                dataj = find(A(:,attr)<=midPoint);
            elseif j == 2
                dataj = find(A(:,attr)>midPoint);
            end
            mod_dataj = length(dataj);
            infoDataj = 0;
            classesj = unique(Y_sorted(dataj));
            for m=1:length(classesj)
                prob_m = sum(ismember(Y_sorted(dataj), classesj{m})) / mod_dataj;
                infoDataj = infoDataj - prob_m * log2(prob_m);
            end
            infoAttrData = infoAttrData + ((mod_dataj / nsamples) * infoDataj);
        end
        infoAttrDatas(i) = infoAttrData;
    end
    minInfoAttrData = min(infoAttrDatas);
    bestSplitPoints{attr} = midPoints{infoAttrDatas == min(infoAttrDatas)};
    % Calculating how much would be gained by branching on attribute
    gainAttrs(attr) = infoData - minInfoAttrData;
end
[sortedGains, sortedAttr] = sortrows(gainAttrs,-1);

end
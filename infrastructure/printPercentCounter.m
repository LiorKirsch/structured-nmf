function  printPercentCounter(i, totalLength)
% printPercentCounter(i, totalLength)
% To be use inside a loop - Prints the percentage passed
% Only updates when a full percentage is hit
    if mod(i, round(totalLength/100)) == 0
        printPercent(round(i*100/totalLength))
    end
end

function printPercent(number)
   prevoiusToPrint = sprintf('%d%%',number-1);
   numberOfChars = length(num2str(prevoiusToPrint));
   fprintf(repmat('\b',1,numberOfChars));
   fprintf('%d%%',number);
end
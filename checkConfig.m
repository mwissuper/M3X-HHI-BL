numInc = 0; numDec = 0; numFollow = 0; numPref = 0;
for i = 1:length(trialArray)
    n = trialArray(i);
    if strcmp(config(n,2),'Increase')
        numInc = numInc + 1;
    elseif strcmp(config(n,2),'Decrease')
        numDec = numDec + 1;
    elseif strcmp(config(n,2),'Follow')
        numFollow = numFollow + 1;
    elseif strcmp(config(n,2),'Preferred')
        numPref = numPref + 1;
    end
end
    
    
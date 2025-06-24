load('workspace_gatecombinations_hillcube_networkC','gatecombs','yyy')
load('dataforprediction_networkC','initial','timedata')

kom=[1,1,1];
for i=1:1:length(timedata)
    timecoursedata{i}=timedata{i}{1};
end


initialconditions=initial/500;
for i=1:1:length(timecoursedata)
    
    normalized_data{i}=[timecoursedata{i}(:,1),(timecoursedata{i}(:,2:end))/500];
    
end

for x=1:1:length(gatecombs(:,1))
    
    parfor i=1:1:length(initialconditions(:,1))
        initial=initialconditions(i,:);
        data=normalized_data(i);
        [cost,ticota] = dataforplotting(gatecombs(x,:),data,yyy,initial,kom);
        FVw(i)=cost;
        stvw(i,:)=ticota{1}(end,:);
    end
    
    OBJ_wit{x}=FVw;
    SD_wit{x}=stvw;
end

save('hillcubepredictionerrors_networkC')
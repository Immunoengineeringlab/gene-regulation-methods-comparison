load('workspace_gatecombinations_rpm_networkA','gatecombs','yyy')
load('dataforprediction_networkA','initial','timedata')

kom=[1,1];
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
        [cost,ticota] = predictionerror(gatecombs(x,:),yyy,data,initial,kom);
        FVr(i)=cost;
        stvr(i,:)=ticota{1}(end,:);
    end
    
    OBJ_rpm{x}=FVr;
    SD_rpm{x}=stvr;
end

save('rpmpredictionerrors_networkA')
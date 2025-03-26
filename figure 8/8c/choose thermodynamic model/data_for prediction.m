load('SaveBest.mat')

sa=Var(50,:); %choose Xth row from Var matrix that corresponds to maximum number of steadystates usually Var(end,:) 


ICS=lhsdesign(100,2)*500;
initial=[zeros(100,3),ICS];
parfor i=1:1:100
    timec = trainingdata(sa,network,initial) %same as "interattractordistance" (objectivefunction file) just the output is timecourse data.
    timedata{i}=timec;
end

save('dataforprediction_networkC')
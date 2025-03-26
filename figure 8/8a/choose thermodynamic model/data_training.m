load('SaveBest.mat')
load('ics.mat','ICS')
sa=Var(50,:); %choose Xth row from Var matrix that corresponds to maximum number of steadystates. usually Var(end,:)
[timecoursedata] = trainingdata(sa,network,ICS) %same as "interattractordistance" (objectivefunction file) just the output is timecourse data.

save("trainingdataA")
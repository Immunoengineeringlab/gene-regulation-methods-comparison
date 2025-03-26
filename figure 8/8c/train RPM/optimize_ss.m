clear all
clc

load('trainingdata','timecoursedata','ICS')

for i=1:1:length(timecoursedata)
    
    data{i}=[timecoursedata{i}(:,1),(timecoursedata{i}(:,2:end))/500]; %normalizing the data to maximum expression level 500.
    
end


kom(1:length(data),1:3)=1; % to assign knockouts. here no KO is assigned to any node.


warning('off')


inputcounts=[5 3 5]; %number of input edges for each node in the order [A B C]

A={'two','three','four','five','six','seven'};

for i=1:1:length(A)
    
    B=A{i};
    yyy(i)=load([B,'inputs'],'mats','G'); %yyy is an object loads the hierarchical canalysing logic gates to be used in RPM or Hillcube simulations
end



for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats;
    UB(q)=length(M(:,1));
end

clear M

lb_cp=[0 0 0  2 2 2  0 0 0 ]';     %lower bound of continuous parameters

ub_cp=[1 1 1  10 10 10 1 1 1 ]';   %upper bound of continuous parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb_dp=[1 1 1]';     %lower bound of discrete parameters representing logic gatges for each node
ub_dp=UB';    %upperbound of discrete parameters representing logic gatges for each node
lb=[lb_dp;lb_cp ];
ub=[ub_dp;ub_cp ];
IntCon = [1 2 3];               % indices of parameters that takes only integer values
nonlcon = [];
options = gaoptimset('UseParallel',true,'display','iter'); % set to "false" if parallel evaluations are not desires
noofvariables=length(ub);

initial=ICS/500;




tic
current_iteration=0;
k=0;
while k<200
    
    current_iteration=current_iteration+1
    
    fun=@(sa) trainrpm_ss(sa,yyy,data,initial,kom);
    tic
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    D=toc
    k=k+1;
    gatecombs_rpm_networkC(k,:)=x;
    Fvalue_rpm_networkC(k,:)=fvalu;
    tiempo(k)=D;
    save('workspace_gatecombinations_rpm_networkC')  %save workspace every 5 iterations
    
end
toc  % evaluate time elapsed
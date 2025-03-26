clear all
clc

load('trainingdata','timecoursedata','ICS')

for i=1:1:length(timecoursedata)
    
    data{i}=[timecoursedata{i}(:,1),(timecoursedata{i}(:,2:end))/500]; %normalizing the data to maximum expression level 500.
    
end


kom(1:length(data),1:4)=1; % to assign knockouts. here no KO is assigned to any node.


warning('off')


inputcounts=[5 3 2 2]; %number of input edges for each node in the order [A B C D]

A={'two','three','four','five'};

for i=1:1:length(A)
    
    B=A{i};
    yyy(i)=load([B,'inputs'],'mats','G'); %yyy is an object loads the hierarchical canalysing logic gates to be used in RPM or Hillcube simulations
    
end


for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats;
    UB(q)=length(M(:,1));
end


lbD=[0 0 0 0];
ubD=[1 1 1 1];

lbK=[0 0 0 0 0 0 0 0 0 0 0 0];
ubK=[1 1 1 1 1 1 1 1 1 1 1 1];

lbn=[2 2 2 2 2 2 2 2 2 2 2 2];
ubn=[5 5 5 5 5 5 5 5 5 5 5 5]*2;


    
initial=ICS/500;




lb_cp=[lbD lbK lbn ]';    %lower bound of continuous parameters
ub_cp=[ubD ubK ubn]';   %upperbound of continuous parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb_dp=[1 1 1 1]';     %lower bound of discrete parameters representing logic gatges for each node
ub_dp=UB';   %upperbound of discrete parameters representing logic gatges for each node
lb=[lb_dp;lb_cp ];
ub=[ub_dp;ub_cp ];
IntCon = [1 2 3 4];               % indices of parameters that takes only integer values
nonlcon = [];
options = gaoptimset('UseParallel',true,'display','iter'); % set to "false" if parallel evaluations are not desires
noofvariables=length(ub);



tic
for I=1:1:200
    
%     data=tmcrsmtrx{I};
    I
    
    fun=@(sa) trainhillcube_ss(sa,data,yyy,initial,kom);
    tic
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    D=toc
    
    tiempo(I)=D;
    
    
    gatecombs_hill_networkB(I,:)=x;
    Fvalue_hill_networkB(I,:)=fvalu;
    save('workspace_gatecombinations_hillcube_networkB')
    I
end
toc


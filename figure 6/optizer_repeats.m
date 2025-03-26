function [x,fvalu] = optizer_repeats(P,cons,network)


noi=length(cons(1,:));

data=P;

xcm=(cons/500);

con=todecimal(network<0);
inputcounts=noi;

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







lb_cp=[0 2]';    %lower bound of parameters
ub_cp=[1 10]';   %upperbound of the parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb_dp=[1]';     %lower bound of parameters
ub_dp=UB';   %upperbound of the parameters
lb=[lb_dp;lb_cp ];
ub=[ub_dp;ub_cp ];
IntCon = [1];               % indices of parameters that takes only integer values
nonlcon = [];
options = gaoptimset('UseParallel',true);%,'display','iter'  % set to "false" if parallel evaluations are not desires
noofvariables=length(ub);





fun=@(sa) rpmoptimize(sa,yyy,data,xcm,con,noi);
[x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);

end


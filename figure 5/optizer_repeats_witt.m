function [x,fvalu] = optizer_repeats_witt(P,cons,network)


noi=length(cons(1,:));

data=P;

xcm=cons/500;

con{1}=todecimal(network<0);
inputcounts=noi;% tbet-6 bcl6 -6 blimp1-5

A={'two','three','four','five'};

for i=1:1:length(A)
    
    B=A{i};
    yyy(i)=load([B,'inputs'],'mats','G'); %yyy is an object loads the hierarchical calalysing logic gates to be used in RPM or Hillcube simulations
    
end


for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats;
    UB(q)=length(M(:,1));
end


lbK=zeros(1,noi);
ubK=ones(1,noi);

lbn=ones(1,noi)+1;
ubn=ones(1,noi)+9;





lb_cp=[lbK lbn]';    %lower bound of parameters
ub_cp=[ubK ubn]';   %upperbound of the parameters
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
options = gaoptimset('UseParallel',true); % ,'display','iter'  % set to "false" if parallel evaluations are not desires
noofvariables=length(ub);





    fun=@(sa) wittmanoptimize(sa,yyy,data,xcm,con,noi);

    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);


end


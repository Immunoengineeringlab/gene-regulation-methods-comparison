clear all
clc



G=[1 1 1 0 0;0 0 0 1 1]; %specify gate -(unit entries in the same row represents conjunction combination and multiple rows represent disjunction combination. for example: G=[1 0;0 1] represents an OR gate

conj=[1 0 0 1 0]; %specify gate negating inputs

%==========================================================================

noi=length(G(:,1)); %number of inputs for RPMgates (number of conjunction terms)
con=G.*conj;



%======================================================================
noiH=length(G(1,:));
decimalVector= 0:1:(2^noiH)-1;
bm= de2bi(decimalVector);

f=[];
for i=1:1:length(bm(:,1))
    q=bm(i,:);
    f(i)=(max(prod(abs((q.^G)-con),2)));
end

F=f;
BM=bm;

%======================================================================

DS=rand(1000,noiH); %randomly sampling the domain of the PBF


for gt=1:1:1000
    
    n=2+3*rand(1,noiH);
    k=rand(1,noiH);
    
    
    for i=1:1:length(DS)
        xc=DS(i,:);
        z(i) = hillcube1(n,k,xc,BM,F);
    end
    
    parametersfordata{gt}=[n;k];
    dataa{gt}=z;
    
end

% surf(x,y,z)
% set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);
% caxis([0 1])
% xlabel('X1')
% ylabel('X2')
% zlabel('BH')

A = [];
b = [];
Aeq = [];
beq = [];

lb=[1 0]';
ub=[10 1]';
IntCon = [];               % indices of parameters that takes only integer values
nonlcon = [];
options = gaoptimset('UseParallel',true); % set to "false" if parallel evaluations are not desiresgaoptimset("Display","iter")
noofvariables=length(ub);

1

tic
for I=1:1:length(dataa)
    
    data=dataa{I};
    
    
    fun=@(sa) RPMgate(sa,data,DS,G,con,noi);
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    
    
    gatecombs(I,:)=x;
    Fvalue(I,:)=fvalu;
    save('workspace_5inputsdata')
    I
end
toc


%to plot histogram in figure 4c

Histogram(Fvalue/1000) % divided by thousand since averaging for 1000 domain sample points (DS)


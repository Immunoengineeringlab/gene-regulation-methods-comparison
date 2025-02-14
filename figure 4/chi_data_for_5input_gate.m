clear all
clc



G{1}=[1 1 1 0 0;0 0 0 1 1];
con{1}=[1 0 0 0 0;0 0 0 1 0];
noi(1)=length(G{1}(1,:));



for T=1:1:length(G)
    g=G{T};
    cono=con{T};
    noio=noi(T);
    %======================================================================
    
    decimalVector= 0:1:(2^noio)-1;
    bm= de2bi(decimalVector);
    conmat=cono.*g;
    f=[];
    for i=1:1:length(bm(:,1))
        q=bm(i,:);
        f(i)=(max(prod(abs((q.^g)-conmat),2)));
    end
    
    F{T}=f;
    BM{T}=bm;
    
end

DS=rand(1000,5); %randomly sampling the domain of the PBF


for gt=1:1:1000
    
    n=2+[3*rand 3*rand 3*rand 3*rand 3*rand];
    k=[rand rand rand rand rand];
    
    
    for i=1:1:length(DS)
        xc=DS(gt,:);
        z(i) = hillcube1(n,k,xc,BM{1},F{1});
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
    
    
    fun=@(sa) fiveinputs(sa,data,DS,G{1},con{1},noi);
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    
    
    gatecombs(I,:)=x;
    Fvalue(I,:)=fvalu;
    save('workspace_5inputsdata')
    I
end
toc


%to plot histogram in figure 4c

Histogram(Fvalue/1000) % divided by thousand since averaging for 1000 domain sample points (DS)


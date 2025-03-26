clear all
clc




%defining gates for Hillcube simulation
G{1}=[1 1]; % AND gate of two inputs. (unit entries in the same row represents conjunction combination and multiple rows represent disjunction combination.
con{1}=[0 0]; % denotes NOT gates of the corresponding inputs

noi(1)=length(G{1}(1,:)); %number of inputs


% hillcube simulation to generate ground truth
%==========================================================================
for T=1:1:length(G)
    g=G{T};
    cono=con{T};
    noio=noi(T);
    
    
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

x=0:0.02:1; %input one 
y=0:0.02:1; %input two


for gt=1:1:1000
    
    n=2+[3*rand 3*rand]; % generates random parameters between two and five  every run
    k=[rand rand]; % generates random parameters between zero and one every run
    
    
    
    
    for i=1:1:length(x)
        for j=1:1:length(y)
            xc=[x(j) y(i)];
            z(i,j) = hillcube1(n,k,xc,BM{1},F{1});
        end
    end
    parametersfordata{gt}=[n;k];
    dataa{gt}=z;
    
end
%==========================================================================

% setting up genetic algorithm optimization for RPM model to reproduce the Hill cube data generated above

A = [];
b = [];
Aeq = [];
beq = [];

lb=[1 0]'; %[n k] 
ub=[10 1]';
IntCon = [];               % indices of parameters that takes only integer values
nonlcon = [];
options = gaoptimset('UseParallel',true); % set to "false" if parallel evaluations are not desiresgaoptimset("Display","iter")
noofvariables=length(ub);



tic 
for I=1:1:length(dataa)

    data=dataa{I};


    fun=@(sa) Chivalue(sa,data);
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);


    gatecombs(I,:)=x;
    Fvalue(I,:)=fvalu;
    save('workspace_chivalues_x1andx2')
    
    I
end
toc


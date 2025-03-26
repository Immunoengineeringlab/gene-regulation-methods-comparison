clear all
clc

warning off

network=[-1,-1,1,-1,1]; % matrix codes for the type of connections to the node. here it denotes 5 inputs with two inhibiting inputs denoted by -1. single row since we have only one node in the network.
PAR=[];
car=[];
for i=1:1:length(network(:,1))
    vect=network(i,:);
    noofspecies=sum(abs(vect));
    A=sum(1:noofspecies);
    PAR=[PAR,rand(1,A),randi([1,4],[1,noofspecies]),rand(1,noofspecies),rand(1,noofspecies)]; % to count the number of parameters required for the thermodynamic model this generates a random parameter vector
    car=[car,rand(1,A),ones([1,noofspecies])*-1,rand(1,noofspecies),rand(1,noofspecies)]; % generates random parameter vector with -1 in the positions where integer values need to be specified for counting the number of molecules of a given species that could bing the promoter
    
end
intcons=find(car<0);

% setting lower and upper bounds for the parameter space
lb=zeros(1,length(PAR));
ub=ones(1,length(PAR));

lb(intcons)=1;
ub(intcons)=4;


A = [];
b = [];
Aeq = [];
beq = [];
IntCon = intcons;   % indices of parameters that takes only integer values (number of binding molecules per species(input) is also a free parameter)
nonlcon = [];
options = gaoptimset('UseParallel',true,"Display","iter",'OutputFcn',@savebestineachgen,'FitnessLimit',0.1,'PopulationSize',10); % set to "false" if parallel evaluations are not desiresgaoptimset("Display","iter")
noofvariables=length(ub);

pd=lhsdesign(1000,length(network))*500; % sampling the concentration levels of each input. 1000 such samples are taken using latin hypercube sampling. arbitrarily high max concentration of 500 is assigned.

tic
for I=1:1:1000
    
    tic
    fun=@(sa) interattractordistance(sa,network,pd);
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    
    
    gatecombs(I,:)=x;
    Fvalue(I,:)=fvalu;
    
    [obj,proba] = PBF_values(x,network,pd); % same function as objective function "interattractordistance". to get the binding probability values at all the sampled input concentrations 
    
    Pbfs(I,:)=proba;
    objs(I,:)=obj;
    
    toc
    save('possible_gates_5inputs')
    I
end
toc



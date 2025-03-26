clear all
clc

warning off

network=[1 -1 -1 1 -1;-1 1 -1 0 0;-1 -1 1 -1 1];  %Network representation nodes in the order [A B C I1 I2]
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
lb=[lb,[0.01 0.01 0.01]];
lb(intcons)=1;
ub=ones(1,length(PAR)+3);
ub(intcons)=4;
lb=[lb,ones(1,3)*0];
ub=[ub,ones(1,3)*5];
A = [];
b = [];
Aeq = [];
beq = [];
IntCon = intcons;   % indices of parameters that takes only integer values (number of binding molecules per species(input) is also a free parameter)
nonlcon = [];
options = gaoptimset('UseParallel',true,"Display","iter",'Generations',700,'OutputFcn',@savebestineachgen,'PopulationSize',50); % set to "false" if parallel evaluations are not desiresgaoptimset("Display","iter")
noofvariables=length(ub);

%setting up initial conditions for inputs I1 and I2

ICS=(dec2bin(0:2^(2)-1)-'0')*10;
initial=(dec2bin(0:2^(2)-1)-'0');
for i=1:1:length(initial(:,1))
    for j=1:1:length(initial(1,:))
        if initial(i,j)==1
            initial(i,j)=450;
        elseif initial(i,j)==0
            initial(i,j)=50;
        end
    end
end

%initial conditions for internal nodes A, B and Cwith stating values of zero

ICS=zeros(4,3);

ICS=[ICS,initial];
% [obj] = interattractordistance(PAR,network,ICS)
save('ics')

tic
for I=1:1:20
    
    
    fun=@(sa) interattractordistance(sa,network,ICS);
    [x,fvalu] = ga(fun,noofvariables,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    
    
    gatecombs(I,:)=x;
    Fvalue(I,:)=fvalu;
    save('workspace_gatecombinations_16ics_4')
    I
end
toc

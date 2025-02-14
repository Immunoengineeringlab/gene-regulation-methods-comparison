function [obj] = interattractordistance(PAR,network,pd)

warning off


% following section is written to convert the linear parameter vector to appropriate matrix form needed for the PBF calculation
% this is necessary sinc optimization algorithm requires vector sampling

x=0;

for i=1:1:length(network(:,1))
    vect=network(i,:);
    noofspecies=sum(abs(vect));
    A=sum(1:noofspecies);
    a=x+1;
    x=x+(A+3*noofspecies);
    par=PAR(a:x);
    [mol_initeractions,q,k,noofbindingmols] = matrixer(vect,par);
    interactionweights{i}=mol_initeractions';
    Qcells{i}=q;
    kcells{i}=k;
    noofbindingmolscell{i}=noofbindingmols;
end

k=1;
intmatrix=interactionweights{k}';
Q=Qcells{k}; % each cell contains the TF-poymerase interaction constants for each TF species
kmx=kcells{k}; % each cell contains the TF-DNA association constants for each TF species
nobms=noofbindingmolscell{k};

noims=length(intmatrix(:,1));

v=1:noims;
tfcombinations=cell(1, (2^noims)-1);
pairs =cell(1, (2^noims)-1);
intweightvect=zeros(1, (2^noims)-1);

%following section computes the possible microstates of the tf
%combinations at the promoter
%polymerase system and computes the interaction weights for each state
k=0;
for i=1:1:noims
    C = nchoosek(v,i);
    for j=1:1:length(C(:,1))
        
        if i==1
            
            k=k+1;
            tfcombinations{k}=C(j,:);
            pairs=[C(j,:),C(j,:)];
            linearindex=sub2ind([noims noims],pairs(:,1),pairs(:,2));
            intweightvect(k)=prod(intmatrix(linearindex));
        else
            
            k=k+1;
            tfcombinations{k}=C(j,:);
            pairs=nchoosek(C(j,:),2);
            linearindex=sub2ind([noims noims],pairs(:,1),pairs(:,2));
            intweightvect(k)=prod(intmatrix(linearindex));
        end
        
    end
end



lp=length(pd(:,1));

proba=zeros(1,lp);


for i=1:1:lp
    
    sa=pd(i,:);
    proba(i)=probabilityofbinding(sa,tfcombinations,intweightvect,Q,kmx,nobms); % calculates binding probability using (1) in He et al.
    
end


% objective is to maximize the spread of probability values in the concentration domain of [0 500] five dimensional space.
%this is important since the choice of some parameter values might lead to probability value of 1 at very low concentrations leading to large
%portion of the domain leading to probability value of 1. Thus to choose
%appropriate parameter values of the thermodynamoc model the objective is
%set up as follows

dist_val=histcounts(proba,'binwidth',0.1); %Probability values are divided into 10 bins and number of concentration samples leading to that bin probability value is counted

s=todecimal(dist_val>50);

obj= 10-sum(s);


end


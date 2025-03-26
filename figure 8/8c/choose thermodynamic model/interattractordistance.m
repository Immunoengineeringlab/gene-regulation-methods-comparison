function [obj] = interattractordistance(PAR,network,ICS)

warning off

% following section is written to convert the linear parameter vector to appropriate matrix form needed for the PBF calculation
% this is necessary since optimization algorithm requires vector sampling

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
    Qcells{i}=q; % each cell contains the TF-poymerase interaction constants for each TF species
    kcells{i}=k; % each cell contains the TF-DNA association constants for each TF species
    noofbindingmolscell{i}=noofbindingmols;
end

%following section computes the possible microstates of the tf
%combinations at the promoter
%polymerase system and computes the interaction weights for each state


for m=1:1:length(network(:,1))
    intmatrix=interactionweights{m};
    noims=length(intmatrix(:,1));
    
    v=1:noims;
    tfcombinations=cell(1, (2^noims)-1);
    intweightvect=zeros(1, (2^noims)-1);
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
    
    odeinputs{1,m}=tfcombinations;
    odeinputs{2,m}=intweightvect;
    
end





initial=ICS;     
concs=[500 500 500];
e=[PAR(end-5) PAR(end-4) PAR(end-3)];
BE=[ PAR(end-2) PAR(end-1) PAR(end)]; %basal expression rate
D=[concs.*e,e]; % calculates  production and degradation rates to be used in the differential equations


gg=0;
for j=1:1:length(initial(:,1))
    
    tstart = tic;
    %length(initial(:,1))
    tspan = [0 6000];
    y0 = initial(j,1:3);
    inputs=initial(j,4:end);
    opts = odeset('Events',@(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs)eventconc(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs,tstart) );
    
    
    [t,y] = ode23s(@toy, tspan, y0,opts,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs);
    
    %timecoursedata{j}=[t,y];
    steadystates(j,:)=y(end,:);
    if t(end)>5500
        gg=1;
    end
    
end


steadystates=steadystates/500;


NOA = uniquesteadystates(steadystates,0.1); % calculates number of different steadystates reached 


obj=     -(intra_CD(steadystates))     - 10*length(NOA(:,1))      +(gg*100000); % maximizes intra cluster distance(intra_CD) with number of unique steadystates(NOA) weighted higher for preference

end


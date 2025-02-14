function [mol_initeractions,q,k,noofbindingmols] = matrixer1(vect,par)


noofspecies=sum(abs(vect));

A=sum(1:noofspecies);

lininteractions=par(1:A);

species_initeractions=ones(noofspecies,noofspecies);
k=0;
for i=1:1:noofspecies
    for j=i:1:noofspecies
        k=k+1;
        if j==i
            species_initeractions(j,i)=1+lininteractions(k)*10; % positive interactions making sure the value is above 1
        else
            species_initeractions(j,i)=lininteractions(k)*1; % negativetive interactions making sure the value is below 1
        end
    end
end

noofbindingmols=par(A+1:A+noofspecies);

mols=[];

for i=1:1:length(noofbindingmols)
    mols=[mols,ones(1,noofbindingmols(i))*i];
end

almo=sum(noofbindingmols);


mol_initeractions=ones(almo,almo);
for i=1:1:almo
    for j=1:1:almo
        if i>j
            mol_initeractions(i,j)=species_initeractions(mols(i),mols(j));
        elseif i==j
            mol_initeractions(i,j)=1;
        end
    end
end

% interactionweights{h}=mol_initeractions;


lininteractions_Q=par(A+1+noofspecies:A+2*noofspecies);

q=[];
dv=vect;
dv(dv==0)=[];
for i=1:1:noofspecies
    qin=lininteractions_Q(i);
%     if dv(i)>0
%         q=[q,1+ones(1,noofbindingmols(i))*qin*2];
%     else
%         q=[q,ones(1,noofbindingmols(i))*qin*0.08];
%     end
    
    if dv(i)>0
        q=[q,1+ones(1,noofbindingmols(i))*(qin*3)]; %for activators- positive interaction between the bound TF and the TRanscription machinary(polymerase) 
    else
        q=[q,ones(1,noofbindingmols(i))*(qin*0.6)]; %for repressors- positive interaction between the bound TF and the TRanscription machinary(polymerase)
    end
    
end
% Qcells{h}=q;

link=par(A+1+2*noofspecies:A+3*noofspecies);
k=[];
% dv=network(h,:);
% dv(dv==0)=[];
for i=1:1:noofspecies
    kin=link(i);
    k=[k,ones(1,noofbindingmols(i))*kin*0.01]; % association constant for TF interaction with the promoter site
end
% kcells{h}=k;
end


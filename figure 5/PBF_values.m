function [obj,proba] = PBF_values(PAR,network,pd)

warning off

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
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};

noims=length(intmatrix(:,1));

v=1:noims;
tfcombinations=cell(1, (2^noims)-1);
pairs =cell(1, (2^noims)-1);
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

TF=zeros(2,(2^noims)-1);


lp=length(pd(:,1));

proba=zeros(1,lp);


for i=1:1:lp

        sa=pd(i,:);
        proba(i)=probabilityofbinding(sa,tfcombinations,intweightvect,Q,kmx,nobms);

end

dist_val=histcounts(proba,'binwidth',0.1);

s=todecimal(dist_val>50);

obj= 10-sum(s);


end


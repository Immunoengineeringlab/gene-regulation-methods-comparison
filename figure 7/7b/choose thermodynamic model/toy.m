function dy = toy(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs)

concs=[y(1),y(2),y(3),y(4),inputs(1),inputs(2)];


k=1;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);
A = D(1)* probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(5)*y(1); %equation 1


k=2;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);
B = D(2)*probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(6)*y(2); %equation 2


k=3;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);
C = D(3)*probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(7)*y(3); %equation 3


k=4;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);
DD = D(4)*probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(8)*y(4); %equation 4


dy = [A,B,C,DD]';


end
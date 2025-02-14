function dy = toy(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs)

concs=[y(1),y(2),inputs(1),inputs(2)];


k=1;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);

A = D(1)* probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(3)*y(1); %equation 1


k=2;
ee=~(network(k,:)==0);
conin=concs(ee);
tfcombinations=odeinputs{1,k};
intweightvect=odeinputs{2,k};
Q=Qcells{k};
kmx=kcells{k};
nobms=noofbindingmolscell{k};
be=BE(k);
B = D(2)*probabilityofbinding(conin,tfcombinations,intweightvect,Q,kmx,nobms,be) - D(4)*y(2); %equation 2






dy = [A,B]';


end
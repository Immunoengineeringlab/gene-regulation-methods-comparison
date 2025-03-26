function [X,isterm,dir] = eventconc(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs,tstart) 




dy = toy(t,y,D,odeinputs,Qcells,kcells,network,noofbindingmolscell,BE,inputs);
dy=(real(dy));
s = abs(dy)<0.00001;


% X(1)=length(dy)-sum(s);
% X(2)=toc(tstart)<2
% 
% 
% isterm = true(size(X));
% dir = zeros(size(X));


X=length(dy)-sum(s);



isterm = 1;
dir = 0;


end
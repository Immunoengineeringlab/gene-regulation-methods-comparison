function [X,isterm,dir] = eventtimecon(t,y,qb,gatecells,LL,kov,I,tstart) 




dy = tfh(t,y,qb,gatecells,LL,kov,I);
dy=(real(dy));
s = abs(dy)<0.00001;


X(1)=length(dy)-sum(s);
X(2)=toc(tstart)<2;

isterm = true(size(X));
dir = zeros(size(X));
end
function [X,isterm,dir] = event(t,y,qb,gatecells,LL,kov,I) 




dy = tfh(t,y,qb,gatecells,LL,kov,I);
s = abs(dy)<0.00001;


X=length(dy)-sum(s);

isterm = 1;
dir = 0;
end
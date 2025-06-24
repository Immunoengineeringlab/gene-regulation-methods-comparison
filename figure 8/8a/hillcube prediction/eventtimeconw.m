function [X,isterm,dir] = eventtimeconw(t,y,D,K,n,BM,F,kov,inputs,tstart)




dy = toyw(t,y,D,K,n,BM,F,kov,inputs);
dy=(real(dy));
s = abs(dy)<0.00001;



X(1)=length(dy)-sum(s);
X(2)=toc(tstart)<2;

isterm = true(size(X));
dir = zeros(size(X));


end
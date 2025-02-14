function [ok] = hillcube1(n,k,xc,bm,f)


% following section calculates the pbf value in equation (17) in the
% manuscript


ee=size(bm);
d=zeros(1,ee(1));

for i=1:1:ee(1)
    
    d1=1;
    for j=1:1:ee(2)
        d1=((bm(i,j)*(((xc(j)^n(j))/(xc(j)^n(j)+k(j)^n(j)))*(1+k(j)^n(j))))+((1-bm(i,j))*(1-(((xc(j)^n(j))/(xc(j)^n(j)+k(j)^n(j)))*(1+k(j)^n(j))))))*d1;
    end
    d(i)=d1;
end

dd=f.*d;
ok=sum(dd);


end


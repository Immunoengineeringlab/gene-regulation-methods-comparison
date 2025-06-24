function [obj] = rpmoptimize(parameters,yyy,data,xcm,con,noi)

gates=parameters(1);

inputcounts=noi;% tbet-6 bcl6 -6 blimp1-5

for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    G=g';
    
end


noih=[size(G,1)];


k=parameters(2);
n=parameters(3);


for i=1:1:length(data)
    
    xc=xcm(i,:);
    r=([(1-xc.^con)+(1-con).*xc]).*G;
    e=r.^G;
    eg=prod(e,2);
    F=(1/noih)*sum(eg);
    z(i)=((F^n)/((F^n)+(k^n)))*(1+k^n);
    
end

obj=sum(abs(data-z));



end
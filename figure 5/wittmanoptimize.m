function [obj] = wittmanoptimize(parameters,yyy,data,xcm,con,noi)


gates=parameters(1);
inputcounts=noi;% tbet-6 bcl6 -6 blimp1-5


K=cell(1);
n=cell(1);



K{1}=parameters(2:2+(noi-1));
n{1}=parameters((2+(noi)):(2+(noi)+(noi-1)));







for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    G{q}=g';
    
end








for T=1:1:length(G)
    g=G{T};
    cono=con{T};
    noio=noi(T);
    %==========================================================================
    
    decimalVector= 0:1:(2^noio)-1;
    bm= de2bi(decimalVector);
    conmat=cono.*g;
    f=[];
    for i=1:1:length(bm(:,1))
        q=bm(i,:);
        f(i)=(max(prod(abs((q.^g)-conmat),2)));
    end
    
    F{T}=f;
    BM{T}=bm;

end


for i=1:1:length(data)
    

    xc=xcm(i,:);

    z(i)=hillcube1(n{1},K{1},xc,BM{1},F{1});
    
end

obj=sum(abs(data-z));



end
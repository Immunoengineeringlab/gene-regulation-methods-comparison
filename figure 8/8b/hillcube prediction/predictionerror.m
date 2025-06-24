function [cost,timecoursedata] = predictionerror(parameters,tmcrsmtrx,yyy,initial,kom)

warning('off')


K=cell(1,4);
n=cell(1,4);

D=parameters(5:8);


K{1}=parameters(9:13);
n{1}=parameters(21:25);


K{2}=parameters(14:16);
n{2}=parameters(26:28);



K{3}=parameters(17:18);
n{3}=parameters(29:30);


K{4}=parameters(19:20);
n{4}=parameters(31:32);


gates=parameters(1:4);
inputcounts=[5 3 2 2];% tbet-6 bcl6 -6 blimp1-5


for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    G{q}=g';
    
end






con{1}=[0 1 0 1 0];
noi(1)=length(G{1}(1,:));


con{2}=[1 0 0];
noi(2)=length(G{2}(1,:));


con{3}=[1 0];
noi(3)=length(G{3}(1,:));

con{4}=[0 1];
noi(4)=length(G{4}(1,:));

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




for j=1:1:length(initial(:,1))
    kov=kom(j,:);
    tspan = [0 600];
    y0 = initial(j,1:4);
    inputs=initial(j,5:end);
    tstart = tic;
    opts = odeset('Events',@(t,y,D,K,n,BM,F,kov,inputs)eventtimeconw(t,y,D,K,n,BM,F,kov,inputs,tstart));
    [t,y] = ode23s(@toyw, tspan, y0,opts,D,K,n,BM,F,kov,inputs);
    timecoursedata{j}=([t,y]);
    
end


for n=1:1:length(initial(:,1))
    
    
    t=timecoursedata{n}(:,1);
    y1=timecoursedata{n}(:,2);
    y2=timecoursedata{n}(:,3);
    y3=timecoursedata{n}(:,4);
    y4=timecoursedata{n}(:,5);

    
    data=tmcrsmtrx{n};
    
    p1=interp1(t,y1,data(:,1));
    p2=interp1(t,y2,data(:,1));
    p3=interp1(t,y3,data(:,1));
    p4=interp1(t,y4,data(:,1));

    %to replace NaN values with steadystate values(NaN shows up in case steadystate is reached before the tpts values)
    p1(isnan(p1))=y1(end);
    p2(isnan(p2))=y2(end);
    p3(isnan(p3))=y3(end);
    p4(isnan(p4))=y4(end);

    % RMSD=sqrt(sum( ((p1-tptvals(1,:))).^2     +   ((p2-tptvals(2,:))).^2   +     ((p3-tptvals(3,:))).^2));
    
    %obj(n)=sum( ((p1-data(:,2))).^2     +   ((p2-data(:,3))).^2 + ((p3-data(:,4))).^2 + ((p4-data(:,5))).^2 ); % chi square estimation
    obj(n)=sum( abs((p1-data(:,2)))     +   abs((p2-data(:,3))) + abs((p3-data(:,4))) + abs((p4-data(:,5))) );
end

cost=sum(obj);

end

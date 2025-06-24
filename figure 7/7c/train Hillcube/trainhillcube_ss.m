function [cost] = trainhillcube_ss(parameters,tmcrsmtrx,yyy,initial,kom)

warning('off')

% following section separates the parameters for each node. Each cell for a parameter corresponds to a node 
K=cell(1,3);
n=cell(1,3);

D=parameters(4:6);


K{1}=parameters(7:11);
n{1}=parameters(20:24);


K{2}=parameters(12:14);
n{2}=parameters(25:27);



K{3}=parameters(15:19);
n{3}=parameters(28:32);

%===========================================================================



gates=parameters(1:3);
inputcounts=[5 3 5];


for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    G{q}=g';
    
end






con{1}=[0 1 1 0 1];
noi(1)=length(G{1}(1,:));


con{2}=[1 0 1];
noi(2)=length(G{2}(1,:));


con{3}=[1 1 0 1 0];
noi(3)=length(G{3}(1,:));



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
    tspan = [0 5000];
    y0 = initial(j,1:3);
    inputs=initial(j,4:end);
    tstart = tic;
    opts = odeset('Events',@(t,y,D,K,n,BM,F,kov,inputs)eventtimeconw(t,y,D,K,n,BM,F,kov,inputs,tstart));
    [t,y] = ode23s(@toyw, tspan, y0,opts,D,K,n,BM,F,kov,inputs);
    timecoursedata{j}=real([t,y]);
    
end

for n=1:1:length(initial(:,1))
    
    if timecoursedata{n}(end,1)<tspan(2)
        y=timecoursedata{n}(end,2:end);
    else
        y=ones(1,3)*2000; %penalty for parameters leading to no steadystates within the simulation time.
    end
    
    data=tmcrsmtrx{n}(end,2:end); % getting the steadystates from the training data
    
    obj(n)=norm(y-data); % euclidian distance between the training data steadystate and the RPM model steadystate 
end

cost=real(sum(obj));

end

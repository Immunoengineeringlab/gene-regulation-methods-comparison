function [cost] = trainrpm_ss(parameters,yyy,tmcrsmtrx,initial,kom)

warning('off')
gates=parameters(1:3);
qb=parameters(4:end);
inputcounts=[5 3 5];

for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    gatecells{q}=g';
    
end


LL=[size(gatecells{1},1) size(gatecells{2},1) size(gatecells{3},1)]; %counts L for each node refer to equation (2) in the manuscript. L is the number of conjunction terms



for j=1:1:length(initial(:,1))
    kov=kom(j,:);
    tspan = [0 5000];
    y0 = initial(j,1:3);
    inputs=initial(j,4:end);
    t=0;
    [xs,~,~] = event(t,y0,qb,gatecells,LL,kov,inputs);
    
    if xs~=0
        tstart = tic;
        opts = odeset('Events',@(t,y,qb,gatecells,LL,kov,inputs)eventtimecon(t,y,qb,gatecells,LL,kov,inputs,tstart),'NonNegative',1);
        [t,y] = ode23s(@tfh, tspan, y0,opts,qb,gatecells,LL,kov,inputs);
        timecoursedata{j}=[t,y];
        
    else
        
        Q=size(tmcrsmtrx{j});
        yd=repmat(y0,Q(1),1);
        timecoursedata{j}=[tmcrsmtrx{j}(:,1),yd];
    end
    
    
end



for n=1:1:length(initial(:,1))
    
    if timecoursedata{n}(end,1)<tspan(2)
        y=real(timecoursedata{n}(end,2:end));
    else
        y=ones(1,3)*2000; %penalty for parameters leading to no steadystates within the simulation time.
    end
    
    data=tmcrsmtrx{n}(end,2:end); % getting the steadystates from the training data
    
    obj(n)=norm(y-data); % euclidian distance between the training data steadystate and the RPM model steadystate 
end
cost=sum(obj);



end
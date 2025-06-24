function [cost,timecoursedata] = predictionerror(parameters,yyy,tmcrsmtrx,initial,kom)

warning('off')
gates=parameters(1:3);
qb=parameters(4:end);
inputcounts=[5 3 5];% tbet-6 bcl6 -6 blimp1-5

for q=1:1:length(inputcounts)
    
    B=inputcounts(q)-1;
    M=yyy(B).mats(gates(q),:);
    w = M(M~=0);
    g=yyy(B).G(:,w);
    gatecells{q}=g';
    
end


LL=[size(gatecells{1},1) size(gatecells{2},1) size(gatecells{3},1)];

%move out
% for i=1:1:length(tmcrsmtrx)
%
%     initial(i,:)=tmcrsmtrx{i}(1,2:end);
%
% end


% gg=0;

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
    
    
    t=timecoursedata{n}(:,1);
    y1=timecoursedata{n}(:,2);
    y2=timecoursedata{n}(:,3);
    y3=timecoursedata{n}(:,4);

    
    data=tmcrsmtrx{n};
    
    p1=interp1(t,y1,data(:,1));
    p2=interp1(t,y2,data(:,1));
    p3=interp1(t,y3,data(:,1));

    %to replace NaN values with steadystate values(NaN shows up in case steadystate is reached before the tpts values)
    p1(isnan(p1))=y1(end);
    p2(isnan(p2))=y2(end);
    p3(isnan(p3))=y3(end);

    % RMSD=sqrt(sum( ((p1-tptvals(1,:))).^2     +   ((p2-tptvals(2,:))).^2   +     ((p3-tptvals(3,:))).^2));
    
    %obj(n)=sum( ((p1-data(:,2))).^2     +   ((p2-data(:,3))).^2 + ((p3-data(:,4))).^2 + ((p4-data(:,5))).^2 ); % chi square estimation
    obj(n)=sum( abs((p1-data(:,2)))     +   abs((p2-data(:,3))) + abs((p3-data(:,4))) );
end

cost=sum(obj);
end
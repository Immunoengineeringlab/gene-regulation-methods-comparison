load('rpmpredictionerrors_networkC')
load('hillcubepredictionerrors_networkC')

for i=1:1:length(OBJ_rpm)
    for j=1:1:length(OBJ_rpm{i})
        tco_rpm(i,j)=OBJ_rpm{i}(1,j);
    end
    
end

for i=1:1:length(OBJ_wit)
    for j=1:1:length(OBJ_wit{i})
        tco_wit(i,j)=OBJ_wit{i}(1,j);
    end
    
end

for i=1:1:length(tco_rpm(1,:))
    X=tco_rpm(:,i);
    X(isnan(X)) = [];
    
    meen_rpm(i)=mean(X);
end

for i=1:1:length(tco_wit(1,:))
    X=tco_wit(:,i);
    X(isnan(X)) = [];
    
    meen_wit(i)=mean(X);
end


Y={meen_rpm',meen_wit'}
rpmc=[240 225 60]/256;
witc=[86 180 233]/256;
colors=[rpmc;witc]
violin(Y,'x',[0 1],'facecolor',colors,'facealpha',1,'edgecolor',[])
xlim([-1 4])

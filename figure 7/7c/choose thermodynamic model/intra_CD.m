function [averagedistance] = intra_CD(Q)

% calculatess the average distance between all pairs of steadystates

X=0;
k=0;
for i=1:1:length(Q(:,1))
    for j=i+1:1:length(Q(:,1))
        k=k+1;
        X=X+norm(Q(i,:)-Q(j,:));
    end
end
averagedistance=X/k;

end


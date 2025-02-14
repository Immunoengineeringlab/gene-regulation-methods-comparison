function [uniquestates] = uniquesteadystates(Q,absdifference)

% calculates unique steadystates from a list of steadystates.

uniquestates=Q;
k=0;
while length(uniquestates(:,1))>k
    k=k+1;
    j=[];
    for i=k+1:1:length(uniquestates(:,1))
        X=sum([abs(uniquestates(k,:)-uniquestates(i,:))>absdifference]);
        if X==0
            j=[j,i];
        end
    end
    j=j';
uniquestates(j,:)=[];
end



end
function [obj] = RPMgate(parameters,data,x,G,con,noi)

n=parameters(1);
k=parameters(2);


for i=1:1:length(x(:,1))

    xcm=x(i,:);
    r=([(1-xcm.^con)+(1-con).*xcm]).*G;
    e=r.^G;
    eg=prod(e,2);
    F=(1/noi)*sum(eg);
    z(i)=((F^n)/((F^n)+(k^n)))*(1+k^n);
    
end

obj=sum(abs(data-z));



end
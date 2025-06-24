function [obj] = Chivalue(parameters,data)

x=0:0.02:1;
y=0:0.02:1;

n=parameters(1);
k=parameters(2);

%RPM simulation

for i=1:1:length(x)
    for j=1:1:length(y)
        a=x(j);
        b=y(i);
        F=(a*b); %PBF
        z(i,j)=((F^n)/((F^n)+(k^n)))*(1+k^n);
    end
end

obj=sum(abs(data-z),"all");

end
close all
% [CII,CI]=min(Fvalue);
%[CII,CI]= min(abs(Fvalue-median(Fvalue))) %to get the location of median Fvalue
[A,B]=sort(Fvalue)
CI=B(72);
%72nd position out of 1000 for min 3B figure since the optimization is
%stochastic user may not be able to generate same image as given in the
%paper but qualitatively similar image will be generated close to lower
%Fvalue. similarly for generating 3D user can choose value close to
%1000th position.
Fvalue(CI)
x=0:0.02:1;
y=0:0.02:1;
hSurface=surf(x,y,dataa{CI})
set(hSurface, 'FaceColor',[86 180 233]/256, 'FaceAlpha',0.5, 'EdgeAlpha', 0);
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);
caxis([0 1])
% xlabel('X1')
% ylabel('X2')
% zlabel('BH')

x=0:0.02:1;
y=0:0.02:1;
n=gatecombs(CI,1);
k=gatecombs(CI,2);

for i=1:1:length(x)
    for j=1:1:length(y)
        a=x(j);
        b=y(i);
        F=(a*b);
        z(i,j)=((F^n)/((F^n)+(k^n)))*(1+k^n);
    end
end

hold on

hSurface=surf(x,y,z) 
set(hSurface, 'FaceColor',[240 225 60]/256, 'FaceAlpha',0.8, 'EdgeAlpha', 0);
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);
caxis([0 1])
% xlabel('X1')
% ylabel('X2')
% zlabel('BH')
% title('1')
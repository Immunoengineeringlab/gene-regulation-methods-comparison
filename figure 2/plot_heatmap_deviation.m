close all
% CI=randi(1000);961
% CI=543
[A,B]=sort(Fvalue)
CI=B(72); %72 is arbitrarily chosen to be a random position close to the minimum value of sorted Fvalue.
%72nd position out of 1000 for min 3C figure since the optimization is
%stochastic user may not be able to generate same image as given in the
%paper but qualitatively similar image will be generated close to lower
%Fvalue. similarly for generating 3E user can choose value close to
%1000th position.

Fvalue(CI)
x=0:0.02:1;
y=0:0.02:1;
hSurface=surf(x,y,dataa{CI})
set(hSurface, 'FaceColor',[0 0.447 0.741], 'FaceAlpha',0.5, 'EdgeAlpha', 0);
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1]);
caxis([0 1])
xlabel('X1')
ylabel('X2')
zlabel('BH')

x=0:0.02:1;
y=0:0.02:1;
n=gatecombs(CI,1);
k=gatecombs(CI,2);

for i=1:1:length(x)
    for j=1:1:length(y)
        a=x(j);
        b=y(i);
        F=(a+b)/2;
        z(i,j)=((F^n)/((F^n)+(k^n)))*(1+k^n);
    end
end

hold on

hSurface=surf(x,y,z) 
set(hSurface, 'FaceColor',[1 1 0], 'FaceAlpha',0.5, 'EdgeAlpha', 0);
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1],'fontsize',16);
caxis([0 1])
xlabel('X1')
ylabel('X2')
zlabel('BH')


rt=abs(z-dataa{CI});
figure


hSurface=surf(x,y,rt)
view(2)
% set(hSurface, 'FaceColor',[0 0.447 0.741], 'FaceAlpha',0.5, 'EdgeAlpha', 0);
set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1],'fontsize',16);
caxis([0 1])
xlabel('X1')
ylabel('X2')
zlabel('BH')
title('1')
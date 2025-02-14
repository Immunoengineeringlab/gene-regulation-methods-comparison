close all
% CI=randi(1000);961
CI=543

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
        F=a*b;
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
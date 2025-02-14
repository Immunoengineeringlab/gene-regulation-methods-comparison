load('workspace_gatecombinations_rpm_networkA','Fvalue_rpm_networkA')

load('workspace_gatecombinations_hillcube_networkA','Fvalue_hill_networkA')

figure(1)
subplot(2,1,1)
histogram(Fvalue_rpm_networA,'FaceColor',[240 225 60]/256,'Normalization','probability')
xlim([0 1])
ylim([0 0.5])
subplot(2,1,2)
histogram(Fvalue_hill_networA,'FaceColor',[86 180 233]/256,'Normalization','probability') 
xlim([0 1])
ylim([0 0.5])
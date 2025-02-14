load('workspace_gatecombinations_rpm_networkB','Fvalue_rpm_networkB')

load('workspace_gatecombinations_hillcube_networkB','Fvalue_hill_networkB')

figure(1)
subplot(2,1,1)
histogram(Fvalue_rpm_networB,'FaceColor',[240 225 60]/256,'Normalization','probability')
xlim([0 1])
ylim([0 0.5])
subplot(2,1,2)
histogram(Fvalue_hill_networB,'FaceColor',[86 180 233]/256,'Normalization','probability') 
xlim([0 1])
ylim([0 0.5])
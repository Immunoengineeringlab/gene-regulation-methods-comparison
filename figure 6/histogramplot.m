fvrpm=load('masterdata_rpm_3inputs.mat','FVALUES')
fvrpm_3 = (fvrpm.FVALUES)/1000
fvhill=load('masterdata_witt_3inputs.mat','FVALUES')
fvhill_3 = (fvhill.FVALUES)/1000


subplot(2,1,1)
histogram(fvrpm_3,'FaceColor',[240 225 60]/256)
xlim
subplot(2,1,2)
histogram(fvhill_3,'FaceColor',[86 180 233]/256) 

fvrpm=load('masterdata_rpm_5inputs.mat','FVALUES')
fvrpm_5 = (fvrpm.FVALUES)/1000
fvhill=load('masterdata_witt_5inputs.mat','FVALUES')
fvhill_5 = (fvhill.FVALUES)/1000

figure
subplot(2,1,1)
histogram(fvrpm_5,'FaceColor',[240 225 60]/256)
xlim
subplot(2,1,2)
histogram(fvhill_5,'FaceColor',[86 180 233]/256)
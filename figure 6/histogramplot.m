fvrpm=load('masterdata_rpm_3inputs.mat','FVALUES')
fvrpm_3 = fvrpm.FVALUES
fvhill=load('masterdata_witt_3inputs.mat','FVALUES')
fvhill_3 = fvhill.FVALUES


subplot(2,1,1)
histogram(fvrpm,'FaceColor',[240 225 60]/256)
xlim
subplot(2,1,2)
histogram(fvhill,'FaceColor',[86 180 233]/256) 

fvrpm=load('masterdata_rpm_5inputs.mat','FVALUES')
fvrpm_5 = fvrpm.FVALUES
fvhill=load('masterdata_witt_5inputs.mat','FVALUES')
fvhill_5 = fvhill.FVALUES

figure
subplot(2,1,1)
histogram(fvrpm,'FaceColor',[240 225 60]/256)
xlim
subplot(2,1,2)
histogram(fvhill,'FaceColor',[86 180 233]/256)
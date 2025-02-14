load('possible_gates_3inputs.mat','Pbfs','pd','network')

for i=1:1:length(Pbfs(:,1))

P=Pbfs(i,:);
cons=pd;

[gatecombs,Fvalue] = optizer_repeats(P,cons,network);

GATES(i,:)=gatecombs;
FVALUES(i)=Fvalue;

i
end

save('masterdata_rpm_3inputs')

clear all

load('possible_gates_5inputs.mat','Pbfs','pd','network')

for i=1:1:length(Pbfs(:,1))

P=Pbfs(i,:);
cons=pd;

[gatecombs,Fvalue] = optizer_repeats(P,cons,network)

GATES(i,:)=gatecombs;
FVALUES(i)=Fvalue;

i
end

save('masterdata_rpm_5inputs')

clear all

load('possible_gates_3inputs.mat','Pbfs','pd','network')

for i=1:1:length(Pbfs(:,1))

P=Pbfs(i,:);
cons=pd;

[gatecombs,Fvalue] = optizer_repeats_witt(P,cons,network);

GATES(i,:)=gatecombs;
FVALUES(i)=Fvalue;

i
end

save('masterdata_witt_3inputs')

clear all

load('possible_gates_5inputs.mat','Pbfs','pd','network')

for i=1:1:length(Pbfs(:,1))

P=Pbfs(i,:);
cons=pd;

[gatecombs,Fvalue] = optizer_repeats_witt(P,cons,network)

GATES(i,:)=gatecombs;
FVALUES(i)=Fvalue;

i
end

save('masterdata_witt_5inputs')
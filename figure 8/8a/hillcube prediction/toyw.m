function dy = toyw(t,y,D,K,n,BM,F,kov,I)

xc=[y(1),y(2),I(1)];

A = D(1)*(   hillcube1(n{1},K{1},xc,BM{1},F{1})      - y(1)); %differential equation for node A


xc=[y(1),y(2),I(2)];

B = D(2)*(    hillcube1(n{2},K{2},xc,BM{2},F{2})     - y(2)); %differential equation for node B






dy = ([A,B].*kov)';


end
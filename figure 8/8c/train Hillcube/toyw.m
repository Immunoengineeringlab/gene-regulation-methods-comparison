function dy = toyw(t,y,D,K,n,BM,F,I)



xc=[y(1),y(2),y(3),I(1),I(2)];

A = D(1)*(   hillcube1(n{1},K{1},xc,BM{1},F{1})      - y(1));  %differential equation for node A






xc=[y(1),y(2),y(3)];

B = D(2)*(    hillcube1(n{2},K{2},xc,BM{2},F{2})     - y(2));  %differential equation for node B







xc=[y(1),y(2),y(3),I(1),I(2)];

C= D(3)*(     hillcube1(n{3},K{3},xc,BM{3},F{3})     - y(3));  %differential equation for node C
 




dy = [A,B,C]';


end
function dy = tfh(t,y,q,g,LL,kov,I)



A = q(5)*(    ((1/(1+(q(1)/((1/LL(1))*(sum(prod([(y(1)),(1-y(2)),(I(1))].^g{1},2)))))^q(3)))*(1+q(1)^q(3)))- y(1));



B =q(6)*(    ((1/(1+(q(2)/((1/LL(2))*(sum(prod([(1-y(1)),(y(2)),(I(2))].^g{2},2)))))^q(4)))*(1+q(2)^q(4)))- y(2));





dy=([A,B].*kov)';


end
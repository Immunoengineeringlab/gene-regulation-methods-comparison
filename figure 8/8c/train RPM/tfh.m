function dy = tfh(t,y,q,g,LL,kov,I)



A = q(7)*(    ((1/(1+(q(1)/((1/LL(1))*(sum(prod([(y(1)),(1-y(2)),(1-y(3)),(I(1)),(1-I(2))].^g{1},2)))))^q(4)))*(1+q(1)^q(4)))- y(1));




B =q(8)*(    ((1/(1+(q(2)/((1/LL(2))*(sum(prod([(1-y(1)),(y(2)),(1-y(3))].^g{2},2)))))^q(5)))*(1+q(2)^q(5)))- y(2));




C = q(9)*(    ((1/(1+(q(3)/((1/LL(3))*(sum(prod([(1-y(1)),(1-y(2)),(y(3)),(1-I(1)),(I(2))].^g{3},2)))))^q(6)))*(1+q(3)^q(6)))- y(3));



dy=([A,B,C].*kov)';

end
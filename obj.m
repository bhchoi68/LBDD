%% Objective function in the Bellman equation

function res = obj(p1,p2)
global beta W0 W1 W2 i1 i2

D0=1-D(p1,p2)-D(p2,p1);
D1=D(p1,p2);
D2=D(p2,p1);
Dsum=D0*W0(i1,i2)+D1*W1(i1,i2)+D2*W2(i1,i2);

res=-(D(p1,p2)*(p1-cost(i1))+beta*Dsum);
end
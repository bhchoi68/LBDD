%%

function res = evol(st,q)

Q1=tranp(q(1));
Q2=tranp(q(2));
dum1=cumsum(Q1(st(1),:));
dum2=cumsum(Q2(st(2),:));
U=rand(2,1);
st1=sum((U(1)>dum1))+1;
st2=sum((U(2)>dum2))+1;

res=[st1 st2];
end
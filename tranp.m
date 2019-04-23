%% L-by-L Markov matrix Pr(w'|w) given q

function res = tranp(q)
global delta L

Del=@(w) 1-(1-delta)^w;
res=zeros(L,L);

if q==0
    res(1,:)=[1 zeros(1,L-1)];
    for i=2:L
        res(i,:)=[zeros(1,i-2) Del(i) 1-Del(i) zeros(1,L-i)];
    end
else
    for i=1:L-1
        res(i,:)=[zeros(1,i-1) Del(i) 1-Del(i) zeros(1,L-i-1)];
    end
    res(L,:)=[zeros(1,L-1) 1];
end
end

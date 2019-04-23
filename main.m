%% Learning-by-doing model, Byunghee Choi

clear; clc;
global L rho kappa l v delta beta W0 W1 W2 i1 i2

%% Parameterization

L=30;
rho=0.85;
kappa=10;
l=15;
v=10;
delta=0.03;
beta=1/1.05;

ini_V=ones(L,L);
ini_p=ones(L,L);
new_V=ini_V;
new_p=ini_p;
dist=1;
tol=1e-4; 
lambda=0.5; % for dampening
num_iter=0;

%% Compute the value and the policy functions of an MPE.

while dist>tol
    W0=ini_V*tranp(0)'*tranp(0)';
    W1=ini_V*tranp(0)'*tranp(1)';
    W2=ini_V*tranp(1)'*tranp(1)';
    for i1=1:L
        for i2=1:L
            new_p(i1,i2)=fminsearch(@(x) obj(x,ini_p(i2,i1)),ini_p(i1,i2));
        end 
    end
    for i1=1:L
        for i2=1:L
            new_V(i1,i2)=-1*obj(new_p(i1,i2),ini_p(i2,i1));
        end
    end
    Vdist=max(max(abs(new_V-ini_V)./(1+abs(new_V))));
    Pdist=max(max(abs(new_p-ini_p)./(1+abs(new_p))));
    dist=max(Vdist,Pdist)
    ini_p=lambda*new_p+(1-lambda)*ini_p;
    ini_V=lambda*new_V+(1-lambda)*ini_V;
    num_iter=num_iter+1
    if num_iter>500
        break
    end
end

%% Starting from state (1,1), Derive the distribution of states after 10, 20, and 30 periods.

time=30; %time-period
num_sim=1000;
st10=zeros(num_sim,2); st20=st10; st30=st10;
dist10=st10; dist20=st20; dist30=st30;
U=rand(time,num_sim);
for sim=1:num_sim
    st=[1,1];%initial state
    sim
    for tt=1:time
        p=[new_p(st(1),st(2)), new_p(st(2),st(1))]; % price choose
        prob=[D(p(1),p(2)),D(p(2),p(1))]; % prob of sale
        q=[(U(tt,sim)<prob(1)), (U(tt,sim)>=prob(1)).*(U(tt,sim)<sum(prob))];
        st=evol(st,q);
        if tt==10
            st10(sim,:)=st;
        elseif tt==20
            st20(sim,:)=st;
        elseif tt==30
            st30(sim,:)=st;
        end
    end
end

for ii=1:L
    for jj=1:L
        dist10(ii,jj)=sum((ii==st10(:,1)).*(jj==st10(:,2)));
        dist20(ii,jj)=sum((ii==st20(:,1)).*(jj==st20(:,2)));
        dist30(ii,jj)=sum((ii==st30(:,1)).*(jj==st30(:,2)));
    end
end
dist10=dist10/num_sim;dist20=dist20/num_sim;dist30=dist30/num_sim;

%% Stationary state distribution.

dist=1;
tol=1e-3;
st=[1,1];
hist_st=[];
tt=0; 
old_dist=dist30;
new_dist=old_dist;
while dist>tol
    tt=tt+1;
    U=rand(1);
    p=[new_p(st(1),st(2)),new_p(st(2),st(1))]; % price choose
    prob=[D(p(1),p(2)),D(p(2),p(1))]; % prob of sale
    q=[(U<prob(1)) (U>=prob(1)).*(U<sum(prob))];
    st=evol(st,q);
    hist_st=[hist_st;st];
    if mod(tt,10000)==0
        check_st=hist_st((tt-7999:tt),:);
        for ii=1:L
            for jj=1:L
                new_dist(ii,jj)=sum((ii==check_st(:,1)).*(jj==check_st(:,2)));
            end
        end
        new_dist=new_dist/10000;
        dist=sum(sum((old_dist-new_dist).^2))
        old_dist=new_dist;
    end
end
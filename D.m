%% Demand: Probability that the firm makes a sale.

function res = D(p1,p2)
global v

res=exp(v-p1)./(1+exp(v-p1)+exp(v-p2));

end
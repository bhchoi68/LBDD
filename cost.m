%% Marginal cost of production (learning curve)

function res = cost(w)
global kappa l rho

eta=log(rho)/log(2);

if w<l
    res=kappa*w^eta;
else
    res=kappa*l^eta;
end
end
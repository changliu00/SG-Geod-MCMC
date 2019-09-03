% function grad = GradESquaredNorm(model, d, t);
function grad = GradESquaredNorm(model, d, t);
% d is in 1:S

T = model.T;
% D = model.D;
V = model.V;
xi = model.xi;

vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0 = sum(vAlphaSD);

aXiSquared = AvK(model.V, model.xi)^2;

grad = 1/(vAlphaSD0*(vAlphaSD0+1)) * ...
       (2*(1-aXiSquared)*model.vMu(1:V, t) + 2*aXiSquared*vAlphaSD(t)*(model.vMu*vAlphaSD));


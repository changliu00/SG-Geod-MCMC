% function derivEsn = DerivESquaredNorm(model, d, j);
function derivEsn = DerivESquaredNorm(model, d, j);
% d here should be in 1:S
T = model.T;
V = model.V;
xi = model.xi;
vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0 = sum(vAlphaSD);

aXiSquared = AvK(V,xi)^2;
vMuTimesVAlpha = model.vMu * vAlphaSD;

esn = ESquaredNorm(model, d);
derivEsn = (1 + 2*(1-aXiSquared)*vAlphaSD(j) + 2*aXiSquared*(vMuTimesVAlpha'*model.vMu(1:V, j)) ...
           - esn*(2*vAlphaSD0+1)) / (vAlphaSD0*(vAlphaSD0 + 1));

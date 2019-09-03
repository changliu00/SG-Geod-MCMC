% function deriv = DerivRhoVAlphaS(model, d, t);
function deriv = DerivRhoVAlphaS(model, d, t);
% d is in 1:S

T = model.T;
xi = model.xi;
V = model.V;
vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0 = sum(vAlphaSD);

aXi = AvK(V,xi);
vMuTimesVAlphaS = model.vMu * vAlphaSD;

esn = ESquaredNorm(model, d);

firstTerm = (model.vMu(1:V, t)/vAlphaSD0 - vMuTimesVAlphaS/vAlphaSD0^2)/sqrt(esn);
secondTerm = -vMuTimesVAlphaS/vAlphaSD0 / (2*esn^(3/2)) * DerivESquaredNorm(model, d, t);

deriv = aXi * dot((firstTerm + secondTerm), model.v(:, model.batchIds(d)));

% function deriv = DerivESquaredNormXi(model,);
function deriv = DerivESquaredNormXi(model);
% size(deriv) is [1 S]
xi = model.xi;
V = model.V;
aXi = AvK(V, xi);
derivAXi = DerivAvK(V, xi);

vAlphaS0s = sum(model.vAlphaS, 1);
sumVAlphasSquared = sum(model.vAlphaS.^2, 1);
%foo = diag(model.vAlphaS' * (model.vMu' * model.vMu) * model.vAlphaS)';
vMuVAlphaVMuVAlpha = sum((model.vAlphaS' * (model.vMu' * model.vMu))' .* model.vAlphaS, 1);
deriv = 2*aXi*derivAXi*(vMuVAlphaVMuVAlpha - sumVAlphasSquared) ./ (vAlphaS0s .* (vAlphaS0s + 1));

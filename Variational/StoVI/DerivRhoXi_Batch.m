% function deriv = DerivRhoXi_Batch(model);
function deriv = DerivRhoXi_Batch(model);

V = model.V;
xi = model.xi;

aXi = AvK(V, xi);
derivAXi = DerivAvK(V, xi);

esns = ESquaredNorm_Batch(model);

vMuTimesVAlphaS = model.vMu * model.vAlphaS;
vAlphaS0s = sum(model.vAlphaS, 1);

deriv = zeros(1, model.S);
derivsOfESquaredNormXi = zeros(1, model.S);
for d = 1:model.S
  derivsOfESquaredNormXi(d) = DerivESquaredNormXi(model, d);
  
  deriv(d) = derivAxi * dot(vMuTimesVAlphaS(:, d), model.v(:, model.batchIds(d))) / vAlphaS0s(d) / sqrt(esns(d)) ...
        - aXi * dot(vMuTimesVAlphaS(:, d), model.v(:, model.batchIds(d))) / 2 / vAlphaS0s(d) / esns(d)^(3/2) * derivsOfESquaredNormXi(d);
end


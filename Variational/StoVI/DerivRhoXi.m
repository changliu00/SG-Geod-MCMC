% function deriv = DerivRhoXi(model);
function deriv = DerivRhoXi(model);

V = model.V;
xi = model.xi;

aXi = AvK(V, xi);
derivAXi = DerivAvK(V, xi);

vAlphaS0s = sum(model.vAlphaS, 1);
esns = ESquaredNorm_Batch(model);
derivESquaredNormXis = DerivESquaredNormXi(model);

vMuTimesVAlphaSDotDoc = sum(model.vAlphaS .* (model.vMu' * model.v(:, model.batchIds)), 1);

%for d = 1:model.S
%  derivESquaredNormXis(d) = DerivESquaredNormXi(model, d);
%end

%deriv = derivAXi * dot(vMuTimesVAlphaS, vd) / vAlphaS0 / sqrt(esn) ...
%       - aXi * dot(vMuTimesVAlphaS, vd) / 2 / vAlphaS0 / esn^(3/2) * derivESquaredNormXi;

deriv = derivAXi * vMuTimesVAlphaSDotDoc ./ (vAlphaS0s .* sqrt(esns)) ...
        - aXi/2 * vMuTimesVAlphaSDotDoc ./ (vAlphaS0s .* esns.^(3/2)) .* derivESquaredNormXis;
      
      

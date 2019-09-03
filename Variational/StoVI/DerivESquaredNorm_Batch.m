% function deriv = DerivESquaredNorm_Batch(model);
function deriv = DerivESquaredNorm_Batch(model);
T = model.T;
V = model.V;
S = model.S;
xi = model.xi;
vAlphaS = model.vAlphaS;
vAlphaS0s = sum(vAlphaS, 1);

aXiSquared = AvK(V,xi)^2;

esns = ESquaredNorm_Batch(model);
vMuTimesVAlphaTimesVMu = model.vAlphaS' * (model.vMu' * model.vMu);  % S by T

% deriv = zeros(T, S);
% for d = 1:S
%   for t = 1:T
%     deriv(t, d) = (1 + 2*(1-aXiSquared)*vAlphaS(t, d) + 2*aXiSquared*vMuTimesVAlphaTimesVMu(d, t) ...
%                    - esns(d)*(2*vAlphaS0s(d)+1)) / (vAlphaS0s(d)*(vAlphaS0s(d) + 1));
%   end
% end

perDocWeights = 1./(vAlphaS0s .* (vAlphaS0s + 1));
deriv = (1 + 2*(1-aXiSquared)*vAlphaS + 2*aXiSquared*vMuTimesVAlphaTimesVMu');
deriv = bsxfun(@minus, deriv,  esns .* (2*vAlphaS0s+1));
deriv = deriv * SparseDiag(perDocWeights);

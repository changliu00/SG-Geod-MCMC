% function esns = ESquaredNorm_Batch(model);
function esns = ESquaredNorm_Batch(model);

T = model.T;
xi = model.xi;
V = model.V;
S = model.S;


aXiSquared = AvK(V,xi)^2;
%vMuVAlphaSVMuVAlphaS = diag(model.vAlphaS' * (model.vMu' * model.vMu) * model.vAlphaS)';
vMuVAlphaSVMuVAlphaS = sum((model.vAlphaS' * (model.vMu' * model.vMu))' .* model.vAlphaS, 1);
vAlphaS0s = sum(model.vAlphaS, 1);
vAlphaSSquares = sum(model.vAlphaS.^2, 1);

% esns = zeros(1, S);
% for d = 1:S
%   vAlphaSD = model.vAlphaS(1:T, d);  % The vAlphaS column for this doc
%   %esns(d) = (vAlphaSD0s(d) + (1-aXiSquared)*vAlphaSSquares(d) ...
%   %          + aXiSquared*vAlphaSD'*vMuTimesVMu*vAlphaSD) / (vAlphaSD0s(d) * (vAlphaSD0s(d)  + 1));
%   esns(d) = (vAlphaSD0s(i) + (1-aXiSquared)*vAlphaSSquares(i) ...
%             + aXiSquared*vAlphaSD'*vMuTimesVMu*vAlphaSD) / (vAlphaSD0s(i) * (vAlphaSD0s(i) + 1));
% end

esns = (vAlphaS0s + (1-aXiSquared)*vAlphaSSquares + aXiSquared*vMuVAlphaSVMuVAlphaS) ./ (vAlphaS0s .* (vAlphaS0s + 1));

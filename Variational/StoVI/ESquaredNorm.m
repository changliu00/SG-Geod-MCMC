% function esn = ESquaredNorm(model, d);
function esn = ESquaredNorm(model, d);
% d is in 1:S

T = model.T;
xi = model.xi;
V = model.V;
vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0 = sum(vAlphaSD);

aXiSquared = AvK(V,xi)^2;
vMuTimesVAlphaS = model.vMu * vAlphaSD;

%topicSquaredNorms = sum(model.vMu.^2, 1); % For Num_Deriv_ correctness
%esn = (vAlphaSD0 + (1-aXiSquared)*sum(vAlphaSD.^2 .* topicSquaredNorms') ...
%        + aXiSquared*(vMuTimesVAlphaS'*vMuTimesVAlphaS)) ...
%        / (vAlphaSD0 * (vAlphaSD0 + 1));
esn = (vAlphaSD0 + (1-aXiSquared)*sum(vAlphaSD.^2) ...
        + aXiSquared*(vMuTimesVAlphaS'*vMuTimesVAlphaS)) ...
        / (vAlphaSD0 * (vAlphaSD0 + 1));

% function deriv = DerivL_alpha(model);
% Returns the gradient of L_alpha w.r.t. alpha.
function deriv = DerivL_alpha(model);
D = model.D;
S = model.S;
kappa1 = model.kappa1;
alpha0 = sum(model.alpha);
vAlphaS0s = sum(model.vAlphaS, 1);  % vAlphaS0 for each document

deriv = sum(psi(model.vAlphaS), 2) - sum(psi(vAlphaS0s)) ...
        + S*psi(alpha0) - S*psi(model.alpha);
deriv = deriv * D / S;


% function deriv = DerivL_vAlphaS(model);
% Returns a matrix containing component-wise derivatives of 
%   d/d{\tilde{\alpha}_d,i} L_{vAlphaS_{d,i}}
function deriv = DerivL_vAlphaS(model);

V = model.V;
S = model.S;
T = model.T;
alpha = model.alpha;
kappa1 = model.kappa1;
alpha0 = sum(alpha);

psiVAlphaS = psi(model.vAlphaS);
derivPsiVAlphaS = psi_n(1, model.vAlphaS);
gammaVAlphaS = gamma(model.vAlphaS);

vAlphaS0s = sum(model.vAlphaS, 1);  % Sum of each document's vAlphaS
gammaVAlphaS0s = gamma(vAlphaS0s);
psiVAlphaS0s = psi(vAlphaS0s);
derivPsiVAlphaS0s = psi_n(1, vAlphaS0s);

%derivsOfRho = zeros(T, S);
%for d = 1:S
%  for t = 1:T
%    derivsOfRho(t, d) = DerivRhoVAlphaS(model, d, t);
%  end
%end
derivsOfRho = DerivRhoVAlphaS_Batch(model); % Compute derivs w.r.t. all vAlphaSs
 
%deriv = zeros(T, S);
%for d = 1:S
%  for t = 1:T
%    deriv(t, d) = ...
%      (alpha(t)-1)*derivPsiVAlphaS(t, d) ...  
%      - derivPsiVAlphaS0s(d)*(alpha0 - T) ... 
%      + kappa1*derivsOfRho(t, d) ... 
%      - (model.vAlphaS(t, d) - 1)*derivPsiVAlphaS(t, d) - psiVAlphaS(t, d) ...
%      + derivPsiVAlphaS0s(d)*(vAlphaS0s(d) - T) + psiVAlphaS0s(d) ...
%      - psiVAlphaS0s(d) + psiVAlphaS(t, d);
%      %+ psiVAlphaS0s(d) - psiVAlphaS(t, d); % 3/19 fix
%  end
%end

% Deriv is T-by-S
deriv = SparseDiag(alpha - 1) * derivPsiVAlphaS ...  % Scale rows by alpha(t)-1
        + kappa1*derivsOfRho ...
        - (model.vAlphaS - 1).*derivPsiVAlphaS - psiVAlphaS ...
        + psiVAlphaS;

addToEachRow = -derivPsiVAlphaS0s*(alpha0 - T) ...
               + derivPsiVAlphaS0s .* (vAlphaS0s - T); 

deriv = bsxfun(@plus, deriv, addToEachRow);



% This is equivalent, and will probably be much faster once if we can
% replace the repmat calls with bsxfun.
%deriv2 = kappa1 * derivsOfRho ...
%         + repmat(derivPsiVAlphaS0s .* (vAlphaS0s - alpha0), T, 1) ...
%         - derivPsiVAlphaS .* (model.vAlphaS - repmat(model.alpha, 1, S));

% function likelihood = L_alpha(model);
% Evaluates L_[alpha], the portion of the log-likelihood dependent on alpha
function likelihood = L_alpha(model);

D = model.D;
S = model.S;
T = model.T;
V = model.V;

alpha = model.alpha;
alpha0 = sum(alpha);

psiVAlphaS = psi(model.vAlphaS);
psiVAlphaS0s = psi(sum(model.vAlphaS, 1));

likelihood = sum(sum( ((alpha - 1) * ones(1, S)) .* psiVAlphaS )) - (alpha0 - T)*sum(psiVAlphaS0s) + ...
             S*LogGamma(alpha0) - S*sum(LogGamma(alpha));
likelihood = likelihood * D / S;

% The expression above should be equivalent to (and more efficient than)
% this old version:
% -----------------
% Sum over documents and the T elements of each alpha
%mySum = 0;
%for d = 1:D
%  vAlphaSD(1:T) = model.vAlphaS(1:T, d);
%  vAlphaSD0 = sum(vAlphaSD);
%  for i = 1:T
%    mySum = mySum + (alpha(i)-1)*(psi(vAlphaSD(i)) - psi(vAlphaSD0));
%  end
%end
%
%likelihood2 = ...
%  mySum + D*log(gamma(alpha0)) - D*sum(log(gamma(alpha)));
% abs(likelihood - likelihood2)

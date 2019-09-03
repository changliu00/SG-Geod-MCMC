% function like = L_vAlphaS(model);
function like = L_vAlphaS(model);

V = model.V;
global D;
global T;
D = model.D;
S = model.S;
T = model.T;
alpha = model.alpha;
kappa1 = model.kappa1;
alpha0 = sum(alpha);

psiVAlphaS = psi(model.vAlphaS);
derivPsiVAlphaS = psi_n(1, model.vAlphaS);

vAlphaS0s = sum(model.vAlphaS, 1);  % Sum of each document's vAlphaS
psiVAlphaS0s = psi(vAlphaS0s);
derivPsiVAlphaS0s = psi_n(1, vAlphaS0s);

alphaMinusOneMatrix = (alpha - 1) * ones(1, S);

%sumOfRhos = 0;
%for d = 1:S
%  sumOfRhos = sumOfRhos + rho(model, d);
%end
sumOfRhos = sum(Rho_Batch(model));

like = sum(sum(alphaMinusOneMatrix .* psiVAlphaS)) - (alpha0 - T)*sum(psiVAlphaS0s) ...
       + S*LogGamma(alpha0) - S*sum(LogGamma(alpha)) ... % WHY NEED THIS TERM???
       + kappa1 * sumOfRhos ...
       - sum(sum((model.vAlphaS - 1) .* psiVAlphaS)) + sum(psiVAlphaS0s .* (vAlphaS0s - T)) ...
       - sum(LogGamma(vAlphaS0s)) + sum(sum(LogGamma(model.vAlphaS)));
       %+ sum(LogGamma(vAlphaS0s)) - sum(sum(LogGamma(model.vAlphaS)));  % 3/19 fix
like = like * D / S;
return;    


%error('Barf');
%like = 0;    
%for d = 1:D
%  like = like + ...
%    sum((alpha-1) .* psiVAlpha(:, d)) - psiVAlpha0s(d)*sum(alpha-1) ...
%    + log(gamma(alpha0)) - sum(log(gamma(alpha))) ...
%    + kappa1 * rho(model, d) ...
%    - sum((model.vAlpha(:, d) - 1) .* psiVAlpha(:, d)) ...
%    + psiVAlpha0s(d)*sum((model.vAlpha(:, d) - 1)) ...
%    + log(gamma(vAlpha0s(d))) - sum(log(gammaVAlpha(:, d)));
%end
       
% term1Fn = @(t, d)((alpha(t)-1)*psiVAlpha(t, d));
% term2Fn = @(t, d)(-(alpha(t)-1)*psiVAlpha0s(d));
% term3Fn = @(t, d)(log(gamma(alpha0)));
% term4Fn = @(t, d)(-log(gamma(alpha(t))));
% term5Fn = @(t, d)(kappa1*rho(model, d));
% term6Fn = @(t, d)(-(model.vAlpha(t, d) - 1) * psiVAlpha(t, d));
% term7Fn = @(t, d)((model.vAlpha(t, d) - 1) * psiVAlpha0s(d));
% term8Fn = @(t, d)(log(gamma(vAlpha0s(d))));
% term9Fn = @(t, d)(-log(gamma(model.vAlpha(t, d))));
% 
% cfLikeTerms = [ ...
%   sum(sum(alphaMinusOneMatrix .* psiVAlpha)), ...
%   - (alpha0 - T)*sum(psiVAlpha0s), ...
%   + D*log(gamma(alpha0)), ...
%   - D*sum(log(gamma(alpha))), ...
%   + kappa1 * sumOfRhos, ...
%   - sum(sum((model.vAlpha - 1) .* psiVAlpha)), ...
%   + sum(psiVAlpha0s .* (vAlpha0s - T)), ...
%   + sum(log(gammaVAlpha0s)), ...
%   - sum(sum(log(gammaVAlpha))) ];
% 
% likeTerms = [ ... 
%   sumOverTD(term1Fn), sumOverTD(term2Fn), ...
%   sumOverD(term3Fn), sumOverTD(term4Fn), ...
%   sumOverD(term5Fn), ...
%   sumOverTD(term6Fn), sumOverTD(term7Fn), ...
%   sumOverD(term8Fn), sumOverTD(term9Fn) ];
% 
%   
% %cfLikeTerms - likeTerms
% 
% function sum = sumOverTD(fn);
% global T;
% global D;
% sum = 0;
% for t = 1:T
%   for d = 1:D
%     sum = sum + fn(t, d);
%   end
% end
% 
% function sum = sumOverD(fn);
% global D;
% sum = 0;
% t = 0;
% for d = 1:D
%   sum = sum + fn(t, d);
% end
% 

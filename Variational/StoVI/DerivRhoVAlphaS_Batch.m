% function derivs = DerivRhoVAlphaS_Batch(model);
% 
% This function does the same thing as DerivRhoVAlphaS, except it computes the 
% derivatives with respect to all vAlphaSs at once.  This removes a lot of
% redundant computation and is significantly faster than calling DerivRhoVAlphaS
% over and over again.
function derivs = DerivRhoVAlphaS_Batch(model);

T = model.T;
xi = model.xi;
V = model.V;
S = model.S;

%vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0s = sum(model.vAlphaS, 1);

aXi = AvK(V,xi);

% A T-by-S matrix containing the derivative of ESN w.r.t. each vAlphaS_t,d
derivsOfESquaredNorm = DerivESquaredNorm_Batch(model);

% A 1-by-S vector containing the expected squared norms of vMu*vAlphaS(:, d)
esns = ESquaredNorm_Batch(model);

vMuDotVd = model.vMu' * model.v(:, model.batchIds);  % T-by-S

%vMuTimesVAlphaS = model.vMu * model.vAlphaS; % One column for each document
%vMuTimesVAlphaSDotVd = sum(vMuTimesVAlphaS .* model.v(:, model.batchIds), 1); % 1-by-S
% Same thing, but way faster and not requiring V*S memory
vMuTimesVAlphaSDotVd = sum(model.vAlphaS .* vMuDotVd, 1);


% derivs = zeros(T, S);
% for d = 1:S
%   for t = 1:T
%     firstTerm = (vMuDotVd(t, d)/vAlphaSD0s(d) - vMuTimesVAlphaSDotVd(d)/vAlphaSD0s(d)^2) / sqrt(esns(d));
%     secondTerm = -vMuTimesVAlphaSDotVd(d)/vAlphaSD0s(d) / (2*esns(d)^(3/2)) * derivsOfESquaredNorm(t, d);
%     derivs(t, d) = aXi * (firstTerm + secondTerm);
% 
%     %firstTerm = (model.vMu(1:V, t)/vAlphaSD0s(d) - vMuTimesVAlphaS(1:V, d)/vAlphaSD0s(d)^2) / sqrt(esns(d));
%     %secondTerm = -vMuTimesVAlphaS(1:V, d)/vAlphaSD0s(d) / (2*esns(d)^(3/2)) * derivsOfESquaredNorm(t, d);
%     %derivs(t, d) = aXi * dot((firstTerm + secondTerm), model.v(:, model.batchIds(d)));
%   end
% end

% Same thing as the nested loop above, except way faster

derivs = vMuDotVd * SparseDiag(1./vAlphaSD0s);

% Subtract a constant from each column
derivs = bsxfun(@minus, derivs, vMuTimesVAlphaSDotVd ./ vAlphaSD0s.^2);

% Divide each column by sqrt(esns(d))
derivs = derivs * SparseDiag(1./sqrt(esns));

s = vMuTimesVAlphaSDotVd ./ vAlphaSD0s ./ (2*esns.^(3/2));
derivs = aXi * (derivs - derivsOfESquaredNorm * SparseDiag(s));


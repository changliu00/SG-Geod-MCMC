% function grad = GradESquaredNorm(model);
% Returns (V, T, D)-dim matrix
%
% TODO: Make a version that computes the gradient for subset of the
% documents -- this will reduce the memory requirements.
%
% TODO: Make a test that verifies that this function does the same thing as
% GradESquaredNorm.
function grad = GradESquaredNorm_Batch(model);

T = model.T;
% D = model.D;
V = model.V;
xi = model.xi;

aXiSquared = AvK(model.V, model.xi)^2;

%vAlphaS0s = sum(model.vAlphaS, 1);  % One entry per document
%vMuTimesVAlphaS = model.vMu * model.vAlphaS; % One column per document
vMuTimesVAlphaS = model.vMu * model.vAlphaS; % One column per document
firstTerm = 2*(1-aXiSquared)*model.vMu;

% grad = zeros(V, T, D);
grad = zeros(V, T, model.S);
%for d = 1:D
for i = 1:model.S
  vAlphaSD = model.vAlphaS(1:T, i);
  vAlphaSD0 = sum(vAlphaSD);
  %secondTerm = (2*aXiSquared) * vMuTimesVAlphaS(:, d) * vAlphaSD';
  %grad(:, :, d) = 1/(vAlphaSD0*(vAlphaSD0+1)) * (firstTerm + secondTerm);
  secondTerm = (2*aXiSquared) * vMuTimesVAlphaS(:, i) * vAlphaSD';
  grad(:, :, i) = 1/(vAlphaSD0*(vAlphaSD0+1)) * (firstTerm + secondTerm);
end

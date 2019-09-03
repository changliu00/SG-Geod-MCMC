% function grad = GradRhoVMu_Batch(model);
%
% Just like GradRhoVMu, except it computes the gradient with respect to each
% document AND topic.  This is way faster than calling GradRhoVMu over and over 
% again, but uses a ton (2000 lbs.) more memory.  Returns a V-by-T-by-D
% matrix.
%
% TODO: Make a version that computes the gradient for a subset of the
% documents.  The caller can partition the documents and call GradRhoVMU_Batch
% on each subset.  This will be less efficient, but will reduce the memory
% requirements.
%
% TODO: Make a test that verifies that this function does the same thing as
% GradRhoVMu.
function grad = GradRhoVMu_Batch(model);

T = model.T;
% D = model.D;
V = model.V;
xi = model.xi;

aXi = AvK(V, xi);
esns = ESquaredNorm_Batch(model);

%gradEsns = GradESquaredNorm_Batch(model);
gradEsns = GradESquaredNorm_Batch(model);

%for t = 1:T
%  gradEsns(:, t) = GradESquaredNorm(model, d, t);
%end
%grad = aXi/vAlphaSD0 * (vAlphaSD(t)*vd/sqrt(esn) - dot(model.vMu*vAlphaSD, vd)/(2*esn^(3/2))* gradEsn);

vAlphaS0s = sum(model.vAlphaS, 1);
%vMuTimesVAlphaSDotV = sum(model.vAlphaS .* (model.vMu' * model.v), 1);
%grad = zeros(V, T, D);
%for d = 1:D

vMuTimesVAlphaSDotV = sum(model.vAlphaS .* (model.vMu' * model.v(:, model.batchIds)), 1);
grad = zeros(V, T, model.S);
for i = 1:model.S
  vAlphaSD = model.vAlphaS(:, i);
  %grad(:, :, d) = aXi/vAlphaS0s(d) * (1/sqrt(esns(d))*model.v(:, d)*vAlphaSD' - vMuTimesVAlphaSDotV(d)/(2*esns(d)^(3/2))* gradEsns(:, :, d));
  grad(:, :, i) = aXi/vAlphaS0s(i) * (1/sqrt(esns(i))*model.v(:, model.batchIds(i))*vAlphaSD' - vMuTimesVAlphaSDotV(i)/(2*esns(i)^(3/2))* gradEsns(:, :, i));
end

%dot(vMu * vAlphaSD, vd) = vAlphaSD' * vMu' * vd  (1-by-T * (T-by-V * V-by-D)
%                       = element d of sum(vAlphaS .* (vMu' * v), 1) 

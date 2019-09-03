% function grad = GradRhoVMu_BatchT(model, d);
%
% Just like GradRhoVMu, except it computes the gradient with respect to each
% topic for the specified document.  This is significantly more efficient than 
% calling GradRhoVMu for each topic individually.  Returns a V-by-T matrix where
% column t contains the gradient of rho(d) with respect to vMu_t.
function grad = GradRhoVMu_BatchT(model, d);
% d is in 1:S

T = model.T;
% D = model.D;
V = model.V;
xi = model.xi;

vAlphaSD = model.vAlphaS(1:T, d);
vAlphaSD0 = sum(vAlphaSD);
vd = model.v(1:V, model.batchIds(d));
aXi = AvK(V, xi);

esn = ESquaredNorm(model, d);
gradEsns = zeros(V, T);
for t = 1:T
  gradEsns(:, t) = GradESquaredNorm(model, d, t);
end
    
%grad = aXi/vAlphaSD0 * (vAlphaSD(t)*vd/sqrt(esn) - dot(model.vMu*vAlphaSD, vd)/(2*esn^(3/2))* gradEsn);
grad = aXi/vAlphaSD0 * (1/sqrt(esn)*vd*vAlphaSD' - dot(model.vMu*vAlphaSD, vd)/(2*esn^(3/2))* gradEsns);

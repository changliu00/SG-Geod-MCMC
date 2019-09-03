% function grad = GradRhoVMu(model, d, t);
function grad = GradRhoVMu(model, d, t);
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
gradEsn = GradESquaredNorm(model, d, t);
    
grad = aXi/vAlphaSD0 * (vAlphaSD(t)*vd / sqrt(esn) - dot(model.vMu*vAlphaSD, vd)/(2*esn^(3/2)) * gradEsn);

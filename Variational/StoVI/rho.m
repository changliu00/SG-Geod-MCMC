% function result = rho(model, d);
function result = rho(model, d);
% d is in 1:S

vAlphaSD = model.vAlphaS(:, d);
vAlphaSD0 = sum(vAlphaSD);
%result = dot(model.vMu * vAlphaSD, model.v(:, model.batchIds(d))) * (AvK(model.V, model.xi) / vAlphaSD0 / sqrt(ESquaredNorm(model, d)));
result = vAlphaSD' * (model.vMu' * model.v(:, model.batchIds(d))) * (AvK(model.V, model.xi) / vAlphaSD0 / sqrt(ESquaredNorm(model, d)));

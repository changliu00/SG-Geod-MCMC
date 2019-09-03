% function result = Rho_Batch(model);
function result = Rho_Batch(model);

% result = vAlphaSD' * (model.vMu' * model.v(:, model.batchIds(d))) * (AvK(model.V, model.xi) / vAlphaSD0 / sqrt(ESquaredNorm(model, d)));

esns = ESquaredNorm_Batch(model);
vAlphaS0s = sum(model.vAlphaS, 1);

result = sum((model.vAlphaS * SparseDiag(1./vAlphaS0s./sqrt(esns))) .* (model.vMu' * model.v(:, model.batchIds)), 1) * AvK(model.V, model.xi);

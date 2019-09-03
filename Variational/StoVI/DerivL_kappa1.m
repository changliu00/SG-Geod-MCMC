% function deriv = DerivL_kappa1(model);
function deriv = DerivL_kappa1(model);
V = model.V;
D = model.D;
S = model.S;

sumOfRhos = sum(Rho_Batch(model));
deriv = -model.D * AvK(model.V, model.kappa1) + sumOfRhos * D / S;

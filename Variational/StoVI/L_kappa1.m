% function likelihood = L_kappa1(model);
function likelihood = L_kappa1(model);
 
D = model.D;
S = model.S;
V = model.V;
k1 = model.kappa1;

sumOfRhos = sum(Rho_Batch(model)) * D / S;
likelihood = D * LogVmfC(V,k1) + k1*sumOfRhos;

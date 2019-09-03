% function grad = DerivL_vMu(model);
function grad = DerivL_vMu(model);
T = model.T;
V = model.V;
xi = model.xi;
V = model.V;
D = model.D;
S = model.S;
k0 = model.kappa0;
k1 = model.kappa1;

%sumOverDocuments = zeros(V, T);
%for d = 1:S
%  sumOverDocuments = sumOverDocuments + GradRhoVMu_BatchT(model, d);
%end
sumOverDocuments = ComputeSumOverDocs(model);
sumOverDocuments = sumOverDocuments * D / S;

% New: vM is a V-by-1 vector, sumOverDocuments is a V-by-T matrix
% Add the first argument to each column of the second argument
grad = bsxfun(@plus, AvK(V,xi) * AvK(V,k0) * xi * model.vM, k1 * sumOverDocuments);

% Project the gradients into the tangent space at each topic
for t = 1:T
  gradForMuT = grad(1:V, t);
  vMuT = model.vMu(1:V, t);
  grad(1:V, t) = gradForMuT - vMuT*(vMuT'*gradForMuT);
end
return;


function s = ComputeSumOverDocs(model);
T = model.T;
S = model.S;
V = model.V;
xi = model.xi;
aXi = AvK(V, xi);
aXiSquared = AvK(model.V, model.xi)^2;
    
esns = zeros(S, 1);  
for d = 1:S
  esns(d) = ESquaredNorm(model, d);
end
vAlphaS0s = sum(model.vAlphaS, 1)';  % Column vector


% For single d:  aXi/vAlphaSD0 /sqrt(esn) * vd * vAlphaSD' 
firstTerm = model.v(:, model.batchIds) * (model.vAlphaS * SparseDiag(aXi ./ (vAlphaS0s .* sqrt(esns))))';

% Weights per document for everything that was in GradESN; from second term in GradRhoVMu_BatchT
% For a single d: aXi/vAlphaSD0 / (2*esn^(3/2)) * dot(model.vMu*vAlphaSD, vd)
perDocWeights = aXi ./ vAlphaS0s ./ (2*esns.^(3/2)) .* sum(model.vAlphaS .* (model.vMu' * model.v(:, model.batchIds)), 1)';

secondTermDocWeights = 1./(vAlphaS0s.*(vAlphaS0s+1))*(2*(1-aXiSquared));
secondTerm = sum(perDocWeights .* secondTermDocWeights) * model.vMu;

% Last term in GradESN times per-doc factors from GradRhoVMu...
thirdTermDocWeights = 1./(vAlphaS0s.*(vAlphaS0s+1)) * 2*aXiSquared;  % From GradESN
rescaledVAlphas = model.vAlphaS * SparseDiag(sqrt(perDocWeights .* thirdTermDocWeights));
thirdTerm = model.vMu * (rescaledVAlphas*rescaledVAlphas');

s = firstTerm - secondTerm - thirdTerm;
return;

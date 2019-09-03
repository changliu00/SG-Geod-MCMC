% function model = VariationalEM(model[, dict[, filePrefix]]]);
function model = VariationalStoEM_CompWithSmp(model, pfx, beginTime, timeInterv)

% SAVE_MODEL = 1;
CONVERGENCE_THRESH = 1e-50;

% Define numerical update rules
vMuRule = struct('paramName', 'vMu', ...
                 'maxIterations', 100, ...
                 'learningRate', 1e-3, ... 
                 'useGradientScaling', 1, ...
                 'gradientMax', 1/sqrt(model.V)/10, ... 
                 'likelihoodFn', @L_vMu, ...
                 'derivFn', @DerivL_vMu, ...
                 'constraintFn', @(x)(NormalizeColumns(x)));
vAlphaSRule = struct('paramName', 'vAlphaS', ...
                    'maxIterations', 300, ...
                    'learningRate', 1e-2, ...
                    'useGradientScaling', 1, ...
                    'gradientMax', 0.05, ...
                    'likelihoodFn', @L_vAlphaS, ...
                    'derivFn', @DerivL_vAlphaS, ...
                    'constraintFn', @(x)(max(x, 1e-5)));
% xiRule = struct('paramName', 'xi', ...
%                 'maxIterations', 100, ...
%                 'learningRate', 1e-2, ...
%                 'useGradientScaling', 1, ...
%                 'gradientMax', 1, ... 
%                 'likelihoodFn', @L_xi, ...
%                 'derivFn', @DerivL_xi, ...
%                 'constraintFn', @(x)(max(x, 1e-5)));
% alphaRule = struct('paramName', 'alpha', ...
%                    'maxIterations', 100, ...
%                    'learningRate', 1e-2, ...
%                    'useGradientScaling', 1, ...
%                    'gradientMax', 0.1, ... 
%                    'likelihoodFn', @L_alpha, ...
%                    'derivFn', @DerivL_alpha, ...
%                    'constraintFn', @(x)(max(x, 1e-5)));

                    
% Helper functions and structures for testing for convergence
isSmall = @(x)(max(max(x)) < CONVERGENCE_THRESH);
parameterNames = {'vMu', 'vM'};

% Initialize some stuff
oldValues = struct();
converged = struct();
model.iteration = 0;

% prepare to write models
stem = sprintf('%s_K-%-3d_sig-%-7.1e_kp0-%-7.1e_kp1-%-7.1e_aph-%-7.1e_tau0-%-2d_gm-%-4.2f', pfx, model.T, model.xi, model.kappa0, model.kappa1, model.alpha(1), model.tau0, model.gamma);
stem( stem == ' ' ) = '_';
rep = 1;
dirname = sprintf('%s_%d', stem, rep);
while exist(dirname, 'dir')
	rep = rep + 1;
	dirname = sprintf('%s_%d', stem, rep);
end
mkdir(dirname);

localTime = 0;
checkPoint = beginTime;
modelNo = 1;
% Go!
while 1 
  model.iteration = model.iteration + 1;  
  disp(sprintf('(Start of iteration %d)', model.iteration));
  % Save old parameter values
  for i = 1:length(parameterNames)
    name = parameterNames{i};
    oldValues.(name) = model.(name);
  end

  % prepare for stochastic VEM
  fullIndices = 1:model.D;
  for s = 1:model.S
      i = randi([s,model.D]);
      model.batchIds(s) = fullIndices(i);
      fullIndices(i) = fullIndices(s);
  end
  omega = 1/(model.iteration + model.tau0)^model.gamma;

  % Perform an EM iteration
  tic;
  model.vAlphaS = model.vAlphaAll(:, model.batchIds);
  model = NumericalUpdate(model, vAlphaSRule);  %model = VEMUpdateVAlphaS(model);
  model.vAlphaAll(:, model.batchIds) = model.vAlphaS;

  model = NumericalUpdate(model, vMuRule);     %model = VEMUpdateVMu(model);
  model.vMu = (1-omega) * oldValues.vMu + omega * model.vMu;
  model.vMu = NormalizeColumns(model.vMu);

  model = VEMUpdateVM(model);
%   model = VEMUpdateM(model);

%   model = NumericalUpdate(model, xiRule);
%	model.xi = (1-omega) * oldValues.xi + omega * model.xi;
%   model = NumericalUpdate(model, alphaRule);
%	model.alpha = (1-omega) * oldValues.alpha + omega * model.alpha;

  %model = VEMUpdateKappa1(model);
  %model = VEMUpdateXi(model);
  %model = VEMUpdateKappa0(model);
  %model = VEMUpdateAlpha(model);

  localTime = toc;
  model.elapsedTime = model.elapsedTime + localTime;
  disp(sprintf('(End of iteration %d, local elapsed time %.3f, total time %.3f)', model.iteration, localTime, model.elapsedTime));

  if model.elapsedTime > checkPoint
	  fprintf('Writing model to disk...');
	  filename = sprintf('%s/%s_%3d_iter-%-5d_Tm-%-.2e.model', dirname, dirname, modelNo, model.iteration, model.elapsedTime);
	  filename( filename == ' ' ) = '_';
	  WriteModel_CompWithSmp(model, filename);
	  fprintf('Done.\n');
	  checkPoint = checkPoint + timeInterv;
	  modelNo = modelNo + 1;
  end
  
  % Check convergence; loop over the parameters
  allParamsConverged = 1;
  for i = 1:length(parameterNames)  
    name = parameterNames{i};
    converged.(name) = isSmall(model.(name) - oldValues.(name));
    if ~converged.(name)
      allParamsConverged = 0;
    end
  end
  
  if allParamsConverged
    fprintf('convergence!\n');
    break;
  end

  if model.iteration >= model.max_iter
	  fprintf('max_iter reached!\n');
	  break;
  end
end
filename = sprintf('%s/%s_iter-last.model', dirname, dirname);
WriteModel_CompWithSmp(model, filename);
disp('Done');


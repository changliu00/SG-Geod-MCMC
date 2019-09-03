% script for comparison with sampling methods

T = 50;
beginTime = 500;
timeInterv = 1000;
max_iter = 200000;
dataFile = '../../Wikipedia-150K/TRD150000_RED2926581_TFIDF1_code.txt';
addpath('.');
addpath('./tools');
addpath('./matlab-host');

model = struct();
model = LoadData_CompWithSmp(model, dataFile);
model = VEMInit_NoVMF_CompWithSmp_wiki(model, T, max_iter);
model = VariationalEM_CompWithSmp(model, 'wiki150K_', beginTime, timeInterv);


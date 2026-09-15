function A10_finalize_reliability_results_test()
% Focused regression coverage for A10_finalize_reliability_results.
% Regressão focada para A10_finalize_reliability_results.
% Exercises two pf_compare_form_mcs.csv variants:
%   1) wide format with Pf_FORM/Pf_MCS_surrogate columns
%   2) long format with Method/Pf/beta/CoV/Error columns

testWideFormat();
testLongFormat();
testZeroReference();

fprintf('A10_finalize_reliability_results_test: OK\n');
end

function testWideFormat()
rootDir = tempname;
cleanupObj = onCleanup(@() cleanupFolder(rootDir)); %#ok<NASGU>
createBaseInputs(rootDir, 0.10, 0.10, 0.12, 0.115);

Tcompare = table( ...
    string({'FORM (surrogate)'; 'MCS (surrogate)'}), ...
    [0.91; 0.92], ...
    [9.1; 9.2], ...
    [0.91; 0.92], ...
    [0.091; 0.092], ...
    [0.11; NaN], ...
    [1.2265; NaN], ...
    [NaN; 0.12], ...
    [NaN; 0.08], ...
    [NaN; 0.003], ...
    'VariableNames', {'Method','Pf','beta','CoV','Error','Pf_FORM','beta_FORM','Pf_MCS_surrogate','CoV_MCS','Error_MCS'});
writeCompareFile(rootDir, Tcompare);

out = A10_finalize_reliability_results(rootDir); %#ok<NASGU>
Tfinal = readtable(fullfile(rootDir, 'out_incremental', 'final_reliability_summary.csv'));

assert(ismember('CoV_Pf', Tfinal.Properties.VariableNames));
assert(ismember('Error_or_StdError', Tfinal.Properties.VariableNames));
assertMetric(Tfinal, 'RS2 empírico', 'AbsError_vs_Pf_ref', 0.0);
assertMetric(Tfinal, 'RS2 empírico', 'RelErrorPct_vs_Pf_ref', 0.0);
assertTextMetric(Tfinal, 'RS2 empírico', 'Bias_vs_Pf_ref', 'coincidente');
assertNaNMetric(Tfinal, 'Pf_ref', 'AbsError_vs_Pf_ref');
assertNaNMetric(Tfinal, 'Pf_ref', 'RelErrorPct_vs_Pf_ref');
assertTextMetric(Tfinal, 'Pf_ref', 'Bias_vs_Pf_ref', 'referência');
assertMetric(Tfinal, 'FORM (surrogate)', 'Pf', 0.11);
assertMetric(Tfinal, 'FORM (surrogate)', 'beta', 1.2265);
assertMetric(Tfinal, 'FORM (surrogate)', 'AbsError_vs_Pf_ref', 0.01);
assertMetric(Tfinal, 'FORM (surrogate)', 'RelErrorPct_vs_Pf_ref', 10.0);
assertTextMetric(Tfinal, 'FORM (surrogate)', 'Bias_vs_Pf_ref', 'superestima');
assertMetric(Tfinal, 'MCS (surrogate)', 'Pf', 0.12);
assertMetric(Tfinal, 'MCS (surrogate)', 'CoV_Pf', 0.08);
assertMetric(Tfinal, 'MCS (surrogate)', 'Error_or_StdError', 0.003);
assertMetric(Tfinal, 'MCS (surrogate)', 'AbsError_vs_Pf_ref', 0.02);
assertMetric(Tfinal, 'MCS (surrogate)', 'RelErrorPct_vs_Pf_ref', 20.0);
assertTextMetric(Tfinal, 'MCS (surrogate)', 'Bias_vs_Pf_ref', 'superestima');
end

function testLongFormat()
rootDir = tempname;
cleanupObj = onCleanup(@() cleanupFolder(rootDir)); %#ok<NASGU>
createBaseInputs(rootDir, 0.10, 0.10, 0.12, 0.115);

Tcompare = table( ...
    string({'FORM'; 'MCS surrogate'}), ...
    [0.105; 0.118], ...
    [1.2536; 1.1840], ...
    [0.02; 0.07], ...
    [0.001; 0.004], ...
    'VariableNames', {'Method','Pf','beta','CoV','Error'});
writeCompareFile(rootDir, Tcompare);

out = A10_finalize_reliability_results(rootDir); %#ok<NASGU>
Tfinal = readtable(fullfile(rootDir, 'out_incremental', 'final_reliability_summary.csv'));

assertMetric(Tfinal, 'FORM (surrogate)', 'Pf', 0.105);
assertMetric(Tfinal, 'RS2 empírico', 'AbsError_vs_Pf_ref', 0.0);
assertMetric(Tfinal, 'RS2 empírico', 'RelErrorPct_vs_Pf_ref', 0.0);
assertTextMetric(Tfinal, 'RS2 empírico', 'Bias_vs_Pf_ref', 'coincidente');
assertNaNMetric(Tfinal, 'Pf_ref', 'AbsError_vs_Pf_ref');
assertNaNMetric(Tfinal, 'Pf_ref', 'RelErrorPct_vs_Pf_ref');
assertTextMetric(Tfinal, 'Pf_ref', 'Bias_vs_Pf_ref', 'referência');
assertMetric(Tfinal, 'FORM (surrogate)', 'CoV_Pf', 0.02);
assertMetric(Tfinal, 'FORM (surrogate)', 'Error_or_StdError', 0.001);
assertMetric(Tfinal, 'FORM (surrogate)', 'AbsError_vs_Pf_ref', 0.005);
assertMetric(Tfinal, 'FORM (surrogate)', 'RelErrorPct_vs_Pf_ref', 5.0);
assertTextMetric(Tfinal, 'FORM (surrogate)', 'Bias_vs_Pf_ref', 'superestima');
assertMetric(Tfinal, 'MCS (surrogate)', 'Pf', 0.118);
assertMetric(Tfinal, 'MCS (surrogate)', 'beta', 1.1840);
assertMetric(Tfinal, 'MCS (surrogate)', 'CoV_Pf', 0.07);
assertMetric(Tfinal, 'MCS (surrogate)', 'Error_or_StdError', 0.004);
assertMetric(Tfinal, 'MCS (surrogate)', 'AbsError_vs_Pf_ref', 0.018);
assertMetric(Tfinal, 'MCS (surrogate)', 'RelErrorPct_vs_Pf_ref', 18.0);
assertTextMetric(Tfinal, 'MCS (surrogate)', 'Bias_vs_Pf_ref', 'superestima');
end

function testZeroReference()
rootDir = tempname;
cleanupObj = onCleanup(@() cleanupFolder(rootDir)); %#ok<NASGU>
createBaseInputs(rootDir, 0.0, 0.0, 0.12, 0.05);

Tcompare = table( ...
    0.11, 1.2265, 0.12, 0.08, 0.003, ...
    'VariableNames', {'Pf_FORM','beta_FORM','Pf_MCS_surrogate','CoV_MCS','Error_MCS'});
writeCompareFile(rootDir, Tcompare);

out = A10_finalize_reliability_results(rootDir); %#ok<NASGU>
Tfinal = readtable(fullfile(rootDir, 'out_incremental', 'final_reliability_summary.csv'));

assertMetric(Tfinal, 'RS2 empírico', 'RelErrorPct_vs_Pf_ref', 0.0);
assertTextMetric(Tfinal, 'RS2 empírico', 'Bias_vs_Pf_ref', 'coincidente');
assertInfMetric(Tfinal, 'FORM (surrogate)', 'RelErrorPct_vs_Pf_ref');
assertInfMetric(Tfinal, 'MCS (surrogate)', 'RelErrorPct_vs_Pf_ref');
end

function createBaseInputs(rootDir, pfRef, pfRs2, pfPck, pfA9)
a5Dir = fullfile(rootDir, 'outputs_a5');
outDir = fullfile(rootDir, 'out_incremental');
mkdir(a5Dir);
mkdir(outDir);

Ta5 = table( ...
    string({'RS2_empirical'; 'A6_PCK'}), ...
    [pfRs2; pfPck], ...
    [betaFromPfTest(pfRs2); betaFromPfTest(pfPck)], ...
    'VariableNames', {'Method','Pf','beta'});
writetable(Ta5, fullfile(a5Dir, 'A5_pf_comparison_summary.csv'));

fid = fopen(fullfile(a5Dir, 'A5_pf_comparison_report.txt'), 'w');
fprintf(fid, 'stub report\n');
fclose(fid);

Tstage0 = table(pfRef, 'VariableNames', {'Pf_ref'});
writetable(Tstage0, fullfile(outDir, 'summary_stage0.csv'));

Tstage4rep = table( ...
    string("OK"), 3, 0.11, 0.10, string("PCK"), 5, ...
    'VariableNames', {'FinalStatus','N_Iterations','Pf_hat_final','Pf_SS_final','bestMethod','nVars'});
writetable(Tstage4rep, fullfile(outDir, 'stage4_report.csv'));

Tstage4hist = table( ...
    [1; 2; 3], ...
    [0.14; 0.12; 0.11], ...
    [0.13; 0.11; 0.10], ...
    [0.90; 0.94; 0.97], ...
    [0.20; 0.15; 0.10], ...
    'VariableNames', {'AL_iter','Pf_hat','Pf_SS','R2_best','LOO_best'});
writetable(Tstage4hist, fullfile(outDir, 'stage4_al_history.csv'));

outSS = struct('Pf', pfA9, 'beta', betaFromPfTest(pfA9), 'CoV', 0.05); %#ok<NASGU>
save(fullfile(outDir, 'A9_subset_simulation_result.mat'), 'outSS');
end

function writeCompareFile(rootDir, Tcompare)
writetable(Tcompare, fullfile(rootDir, 'out_incremental', 'pf_compare_form_mcs.csv'));
end

function assertMetric(T, methodName, columnName, expectedValue)
idx = strcmp(string(T.Method), string(methodName));
assert(any(idx), 'Método não encontrado: %s', methodName);

actualValue = T.(columnName)(find(idx, 1, 'first'));
assert(abs(actualValue - expectedValue) < 1e-12, ...
    'Valor inesperado para %s / %s: %.15g ~= %.15g', ...
    methodName, columnName, actualValue, expectedValue);
end

function assertTextMetric(T, methodName, columnName, expectedValue)
idx = strcmp(string(T.Method), string(methodName));
assert(any(idx), 'Método não encontrado: %s', methodName);

actualValue = string(T.(columnName)(find(idx, 1, 'first')));
assert(actualValue == string(expectedValue), ...
    'Texto inesperado para %s / %s: %s ~= %s', ...
    methodName, columnName, actualValue, expectedValue);
end

function assertNaNMetric(T, methodName, columnName)
idx = strcmp(string(T.Method), string(methodName));
assert(any(idx), 'Método não encontrado: %s', methodName);

actualValue = T.(columnName)(find(idx, 1, 'first'));
assert(isnan(actualValue), 'Esperado NaN para %s / %s.', methodName, columnName);
end

function assertInfMetric(T, methodName, columnName)
idx = strcmp(string(T.Method), string(methodName));
assert(any(idx), 'Método não encontrado: %s', methodName);

actualValue = T.(columnName)(find(idx, 1, 'first'));
assert(isinf(actualValue), 'Esperado Inf para %s / %s.', methodName, columnName);
end

function beta = betaFromPfTest(pf)
epsv = 1e-15;
pf = min(max(pf, epsv), 1 - epsv);
beta = -norminv(pf);
end

function cleanupFolder(rootDir)
if isfolder(rootDir)
    rmdir(rootDir, 's');
end
end

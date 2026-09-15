function out = A10_finalize_reliability_results(work_dir)
% A10_FINALIZE_RELIABILITY_RESULTS
% Consolida os resultados finais de confiabilidade para uso em dissertação.
%
% Entradas esperadas do workflow:
%   - outputs_a5/A5_pf_comparison_summary.csv
%   - outputs_a5/A5_pf_comparison_report.txt (opcional)
%   - pf_compare_form_mcs.csv
%   - summary_stage0.csv
%   - out_incremental/stage4_report.csv
%   - out_incremental/stage4_al_history.csv
%   - out_incremental/A9_subset_simulation_result.mat
%
% Saídas:
%   - out_incremental/final_reliability_summary.csv
%   - out_incremental/final_reliability_report.txt
%   - out_incremental/final_reliability_comparison.png

if nargin < 1 || isempty(work_dir)
    work_dir = 'C:\Kazuo-Script';
end

%% Paths
out_dir  = fullfile(work_dir, 'out_incremental');
a5_dir   = fullfile(work_dir, 'outputs_a5');

f_a5_summary  = locateRequiredFile('A5_pf_comparison_summary.csv', {a5_dir, out_dir, work_dir});
f_a5_report   = locateOptionalFile('A5_pf_comparison_report.txt', {a5_dir, out_dir, work_dir});
f_pf_compare  = locateRequiredFile('pf_compare_form_mcs.csv', {out_dir, a5_dir, work_dir});
f_stage0      = locateRequiredFile('summary_stage0.csv', {out_dir, a5_dir, work_dir});
f_stage4_rep  = locateRequiredFile('stage4_report.csv', {out_dir, a5_dir, work_dir});
f_stage4_hist = locateRequiredFile('stage4_al_history.csv', {out_dir, a5_dir, work_dir});
f_a9_result   = locateRequiredFile('A9_subset_simulation_result.mat', {out_dir, a5_dir, work_dir});

f_out_csv = fullfile(out_dir, 'final_reliability_summary.csv');
f_out_txt = fullfile(out_dir, 'final_reliability_report.txt');
f_out_png = fullfile(out_dir, 'final_reliability_comparison.png');

dpi = 220;

if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% Read inputs
Ta5        = readtable(f_a5_summary);
Tcompare   = readtable(f_pf_compare);
Tstage0    = readtable(f_stage0);
Tstage4rep = readtable(f_stage4_rep);
Tstage4hist = readtable(f_stage4_hist);
S9         = load(f_a9_result);
outSS      = extractSubsetSimulationOutput(S9);

if ~isempty(f_a5_report) && isfile(f_a5_report)
    txtA5 = fileread(f_a5_report); %#ok<NASGU>
end

%% Extract main reliability metrics
Pf_ref = getStage0PfRef(Tstage0);
beta_ref = pf2beta(Pf_ref);

Pf_rs2 = getMethodMetric(Ta5, {'RS2_empirical','RS2 empirical','RS2'}, {'Pf'});
beta_rs2 = getMethodMetric(Ta5, {'RS2_empirical','RS2 empirical','RS2'}, {'beta','Beta'});
if isnan(beta_rs2)
    beta_rs2 = pf2beta(Pf_rs2);
end

Pf_pck = getMethodMetric(Ta5, {'PCK','A6_PCK','PCK_final'}, {'Pf'});
beta_pck = getMethodMetric(Ta5, {'PCK','A6_PCK','PCK_final'}, {'beta','Beta'});
if isnan(beta_pck)
    beta_pck = pf2beta(Pf_pck);
end

[Pf_form, colPfForm] = extractCompareMetric(Tcompare, {'Pf_FORM','FORM_Pf','PfForm','Pf_FORM_surrogate'}, {'FORM','FORM surrogate','FORM (surrogate)','FORM_surrogate'}, {'Pf'});
[beta_form, colBetaForm] = extractCompareMetric(Tcompare, {'beta_FORM','Beta_FORM','FORM_beta','betaForm','beta_FORM_surrogate'}, {'FORM','FORM surrogate','FORM (surrogate)','FORM_surrogate'}, {'beta','Beta'});
if isnan(beta_form)
    beta_form = pf2beta(Pf_form);
end
[cov_form, colCovForm] = extractCompareMetric(Tcompare, {'CoV_FORM','cov_FORM','CoV_Pf_FORM'}, {'FORM','FORM surrogate','FORM (surrogate)','FORM_surrogate'}, {'CoV','cov'});
[err_form, colErrForm] = extractCompareMetric(Tcompare, {'Error_FORM','StdErr_FORM','SE_FORM','AbsError_FORM'}, {'FORM','FORM surrogate','FORM (surrogate)','FORM_surrogate'}, {'Error','StdErr','SE'});

[Pf_mcs, colPfMcs] = extractCompareMetric( ...
    Tcompare, ...
    {'Pf_MCS','Pf_MCS_surrogate','Pf_surrogate_on_RS2','Pf_surrogate_on_RS2_calibrated','Pf_MC','Pf_MonteCarlo'}, ...
    {'MCS','MCS surrogate','MCS (surrogate)','MCS_surrogate','Monte Carlo','MonteCarlo'}, ...
    {'Pf'});
[beta_mcs, colBetaMcs] = extractCompareMetric(Tcompare, {'beta_MCS','Beta_MCS','MCS_beta','beta_MC'}, {'MCS','MCS surrogate','MCS (surrogate)','MCS_surrogate','Monte Carlo','MonteCarlo'}, {'beta','Beta'});
if isnan(beta_mcs)
    beta_mcs = pf2beta(Pf_mcs);
end
[cov_mcs, colCovMcs] = extractCompareMetric( ...
    Tcompare, ...
    {'CoV_MCS','cov_MCS','CoV_Pf_MCS','CoV_surrogate_on_RS2','cov_surrogate_on_RS2','CoV_Pf_surrogate_on_RS2'}, ...
    {'MCS','MCS surrogate','MCS (surrogate)','MCS_surrogate','Monte Carlo','MonteCarlo'}, ...
    {'CoV','cov'});
[err_mcs, colErrMcs] = extractCompareMetric( ...
    Tcompare, ...
    {'Error_MCS','StdErr_MCS','SE_MCS','AbsError_MCS','Error_surrogate_on_RS2'}, ...
    {'MCS','MCS surrogate','MCS (surrogate)','MCS_surrogate','Monte Carlo','MonteCarlo'}, ...
    {'Error','StdErr','SE'});

Pf_ss = getStructMetric(outSS, {'Pf','pf'});
beta_ss = getStructMetric(outSS, {'beta','Beta'});
if isnan(beta_ss)
    beta_ss = pf2beta(Pf_ss);
end
cov_ss = getStructMetric(outSS, {'CoV','cov','CoV_Pf'});

%% Final comparison table
Method = string({
    'Pf_ref'
    'RS2 empírico'
    'PCK'
    'FORM (surrogate)'
    'MCS (surrogate)'
    'Subset Simulation (A9 surrogate)'
});

Pf = [
    Pf_ref
    Pf_rs2
    Pf_pck
    Pf_form
    Pf_mcs
    Pf_ss
];

Beta = [
    beta_ref
    beta_rs2
    beta_pck
    beta_form
    beta_mcs
    beta_ss
];

CoV = [
    NaN
    NaN
    NaN
    cov_form
    cov_mcs
    cov_ss
];

ErrorOrStd = [
    NaN
    NaN
    NaN
    err_form
    err_mcs
    NaN
];

MetricDetail = string({
    'referência stage 0'
    'estimativa empírica em RS2'
    'metamodelo PCK do fluxo principal'
    buildMetricDetail(colPfForm, colBetaForm, colCovForm, colErrForm)
    buildMetricDetail(colPfMcs, colBetaMcs, colCovMcs, colErrMcs)
    'resultado A9_subset_simulation_result.mat'
});

AbsErrorPfRef = NaN(size(Pf));
RelErrorPctPfRef = NaN(size(Pf));
RelDirection = strings(size(Method));
RelDirection(:) = missing;
RelDirection(1) = "referência";
for i = 2:numel(Method)
    if isnan(Pf(i))
        continue;
    end
    AbsErrorPfRef(i) = abs(Pf(i) - Pf_ref);
    RelErrorPctPfRef(i) = computeRelativeErrorPct(AbsErrorPfRef(i), Pf_ref);
    RelDirection(i) = classifyRelativeBias(Pf(i), Pf_ref);
end

Tfinal = table( ...
    Method, Pf, Beta, CoV, ErrorOrStd, AbsErrorPfRef, RelErrorPctPfRef, RelDirection, MetricDetail, ...
    'VariableNames', {'Method','Pf','beta','CoV_Pf','Error_or_StdError','AbsError_vs_Pf_ref','RelErrorPct_vs_Pf_ref','Bias_vs_Pf_ref','SourceDetail'});

writetable(Tfinal, f_out_csv);

%% Supporting stage 4 diagnostics
finalStatus = string(getFirstTableValue(Tstage4rep, {'FinalStatus','Status'}));
nIterStage4 = getFirstTableValue(Tstage4rep, {'N_Iterations','NIterations','nIter','Iterations'});
Pf_hat_final = getFirstTableValue(Tstage4rep, {'Pf_hat_final','Pf_hat','PfHatFinal'});
Pf_SS_final = getFirstTableValue(Tstage4rep, {'Pf_SS_final','Pf_SS','PfSSFinal'});
bestMethod4 = string(getFirstTableValue(Tstage4rep, {'bestMethod','BestMethod','Method'}));
nVars4 = getFirstTableValue(Tstage4rep, {'nVars','NVars','n_variables'});

[histIter, colHistIter] = getHistoryVector(Tstage4hist, {'AL_iter','Iter','Iteration','ALIteration'});
[histPfHat, colHistPfHat] = getHistoryVector(Tstage4hist, {'Pf_hat','PfHat'});
[histPfSS, colHistPfSS] = getHistoryVector(Tstage4hist, {'Pf_SS','PfSS'});
[histR2, colHistR2] = getHistoryVector(Tstage4hist, {'R2_best','R2','R2Best'});
[histLOO, colHistLOO] = getHistoryVector(Tstage4hist, {'LOO_best','LOO','LOOBest','Q2_LOO'});

%% 4-panel figure
f = figure('Color','w','Position',[100 100 1500 950],'Visible','off');
figureCleanup = onCleanup(@() safeCloseFigure(f)); %#ok<NASGU>
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% Panel 1: final Pf comparison
nexttile;
plotComparisonBars(Method, Pf, [0.24 0.48 0.78]);
ylabel('P_f');
title('Comparação final de P_f');

% Panel 2: final beta comparison
nexttile;
plotComparisonBars(Method, Beta, [0.86 0.42 0.20]);
ylabel('\beta');
title('Comparação final de \beta');

% Panel 3: relative error vs Pf_ref
nexttile;
errMethods = Method;
errVals = RelErrorPctPfRef;
plotComparisonBars(errMethods(2:end), errVals(2:end), [0.35 0.64 0.35]);
ylabel('Erro relativo [%]');
title('Erro relativo em relação a P_f_{ref}');
hold on;
yline(0, '--k', 'LineWidth', 1.0);
hold off;

% Panel 4: stage 4 evolution
nexttile;
hasLeft = false;
hasRight = false;
if ~isempty(histIter)
    if ~isempty(histPfHat)
        [iterPfHat, valPfHat] = buildPlotSeries(histIter, histPfHat);
    else
        iterPfHat = [];
        valPfHat = [];
    end
    if ~isempty(histPfSS)
        [iterPfSS, valPfSS] = buildPlotSeries(histIter, histPfSS);
    else
        iterPfSS = [];
        valPfSS = [];
    end
    if ~isempty(histR2)
        [iterR2, valR2] = buildPlotSeries(histIter, histR2);
    else
        iterR2 = [];
        valR2 = [];
    end
    if ~isempty(histLOO)
        [iterLOO, valLOO] = buildPlotSeries(histIter, histLOO);
    else
        iterLOO = [];
        valLOO = [];
    end

    if ~isempty(valPfHat)
        yyaxis left;
        plot(iterPfHat, valPfHat, '-o', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
        hasLeft = true;
    end
    if ~isempty(valPfSS)
        if ~hasLeft
            yyaxis left;
            hold on;
        end
        plot(iterPfSS, valPfSS, '-s', 'LineWidth', 1.5, 'MarkerSize', 5);
        hasLeft = true;
    end
    if hasLeft
        ylabel('P_f');
    end

    if ~isempty(valR2) || ~isempty(valLOO)
        yyaxis right;
        hold on;
        if ~isempty(valR2)
            plot(iterR2, valR2, '-^', 'LineWidth', 1.5, 'MarkerSize', 5);
        end
        if ~isempty(valLOO)
            plot(iterLOO, valLOO, '-d', 'LineWidth', 1.5, 'MarkerSize', 5);
        end
        ylabel('Qualidade do surrogate');
        hasRight = true;
    end

    if hasLeft || hasRight
        grid on;
        xlabel('Iteração AL');
        title('Diagnósticos de stage 4');
        legendStrings = buildStage4Legend(~isempty(valPfHat), ~isempty(valPfSS), ~isempty(valR2), ~isempty(valLOO));
        legend(legendStrings, 'Location', 'best');
    else
        text(0.1, 0.5, 'stage4\_al\_history.csv sem colunas suportadas', 'FontSize', 11);
        axis off;
    end
else
    text(0.1, 0.5, 'stage4\_al\_history.csv sem coluna de iteração', 'FontSize', 11);
    axis off;
end

annotationText = sprintf('Histórico: iter=%s | Pf_{hat}=%s | Pf_{SS}=%s | R2=%s | LOO=%s', ...
    fallbackLabel(colHistIter), fallbackLabel(colHistPfHat), fallbackLabel(colHistPfSS), fallbackLabel(colHistR2), fallbackLabel(colHistLOO));
sgtitle({'Consolidação final de confiabilidade', annotationText}, 'FontWeight', 'bold');

exportgraphics(f, f_out_png, 'Resolution', dpi);

%% Text report
fid = fopen(f_out_txt, 'w');
if fid <= 0
    error('Nao foi possivel criar o relatorio final: %s', f_out_txt);
end
cleanupObj = onCleanup(@() safeCloseFile(fid)); %#ok<NASGU>
methodWidth = max(32, max(strlength(Tfinal.Method)) + 2);

fprintf(fid, 'RELATORIO FINAL DE CONFIABILIDADE\n');
fprintf(fid, '=================================\n\n');
fprintf(fid, 'Data/hora: %s\n\n', datestr(now));

fprintf(fid, '1) Arquivos de entrada consolidados\n');
fprintf(fid, '   A5 summary      : %s\n', f_a5_summary);
fprintf(fid, '   FORM/MCS compare: %s\n', f_pf_compare);
fprintf(fid, '   Stage 0 summary : %s\n', f_stage0);
fprintf(fid, '   Stage 4 report  : %s\n', f_stage4_rep);
fprintf(fid, '   Stage 4 history : %s\n', f_stage4_hist);
fprintf(fid, '   A9 result       : %s\n\n', f_a9_result);

fprintf(fid, '2) Comparacao final principal (referencia = Pf_ref)\n');
fprintf(fid, '   %-*s %12s %12s %12s %12s %14s %14s %14s\n', methodWidth, 'Metodo', 'Pf', 'beta', 'CoV', 'Err/SE', '|erro abs|', 'erro rel [%]', 'bias');
fprintf(fid, '   %s\n', repmat('-',1, methodWidth + 94));
for i = 1:height(Tfinal)
    fprintf(fid, '   %-*s %12s %12s %12s %12s %14s %14s %14s\n', ...
        methodWidth, char(Tfinal.Method(i)), ...
        formatNumber(Tfinal.Pf(i), '%.6g'), ...
        formatNumber(Tfinal.beta(i), '%.6f'), ...
        formatNumber(Tfinal.CoV_Pf(i), '%.6g'), ...
        formatNumber(Tfinal.Error_or_StdError(i), '%.6g'), ...
        formatNumber(Tfinal.AbsError_vs_Pf_ref(i), '%.6g'), ...
        formatNumber(Tfinal.RelErrorPct_vs_Pf_ref(i), '%.3f'), ...
        char(string(Tfinal.Bias_vs_Pf_ref(i))));
end
fprintf(fid, '\n');

fprintf(fid, '3) Leitura do arquivo pf_compare_form_mcs.csv\n');
fprintf(fid, '   Pf_FORM selecionado a partir de       : %s\n', fallbackLabel(colPfForm));
fprintf(fid, '   beta_FORM selecionado a partir de     : %s\n', fallbackLabel(colBetaForm));
fprintf(fid, '   CoV/erro FORM selecionado a partir de : %s / %s\n', fallbackLabel(colCovForm), fallbackLabel(colErrForm));
fprintf(fid, '   Pf_MCS selecionado a partir de        : %s\n', fallbackLabel(colPfMcs));
fprintf(fid, '   beta_MCS selecionado a partir de      : %s\n', fallbackLabel(colBetaMcs));
fprintf(fid, '   CoV/erro MCS selecionado a partir de  : %s / %s\n\n', fallbackLabel(colCovMcs), fallbackLabel(colErrMcs));

fprintf(fid, '4) Leitura interpretativa frente a Pf_ref\n');
writeInterpretation(fid, 'RS2 empírico', Pf_rs2, Pf_ref);
writeInterpretation(fid, 'PCK', Pf_pck, Pf_ref);
writeInterpretation(fid, 'FORM (surrogate)', Pf_form, Pf_ref);
writeInterpretation(fid, 'MCS (surrogate)', Pf_mcs, Pf_ref);
writeInterpretation(fid, 'Subset Simulation (A9 surrogate)', Pf_ss, Pf_ref);
fprintf(fid, '\n');

fprintf(fid, '5) Diagnosticos de stage 4 (apoio)\n');
fprintf(fid, '   FinalStatus   = %s\n', char(finalStatus));
fprintf(fid, '   N_Iterations  = %s\n', fmtValue(nIterStage4));
fprintf(fid, '   Pf_hat_final  = %s\n', fmtValue(Pf_hat_final));
fprintf(fid, '   Pf_SS_final   = %s\n', fmtValue(Pf_SS_final));
fprintf(fid, '   bestMethod    = %s\n', char(bestMethod4));
fprintf(fid, '   nVars         = %s\n', fmtValue(nVars4));
fprintf(fid, '   Colunas historico usadas: iter=%s | Pf_hat=%s | Pf_SS=%s | R2=%s | LOO=%s\n\n', ...
    fallbackLabel(colHistIter), fallbackLabel(colHistPfHat), fallbackLabel(colHistPfSS), fallbackLabel(colHistR2), fallbackLabel(colHistLOO));

fprintf(fid, '6) Saidas geradas\n');
fprintf(fid, '   %s\n', f_out_csv);
fprintf(fid, '   %s\n', f_out_txt);
fprintf(fid, '   %s\n', f_out_png);

%% Console output
fprintf('\n=== A10: consolidacao final concluida ===\n');
disp(Tfinal);
fprintf('Relatorio : %s\n', f_out_txt);
fprintf('Figura    : %s\n', f_out_png);

%% Output struct
out = struct();
out.csv = f_out_csv;
out.txt = f_out_txt;
out.png = f_out_png;
out.table = Tfinal;
out.Pf_ref = Pf_ref;
out.Pf_rs2 = Pf_rs2;
out.Pf_pck = Pf_pck;
out.Pf_form = Pf_form;
out.Pf_mcs = Pf_mcs;
out.Pf_ss = Pf_ss;
out.stage4_status = finalStatus;
out.stage4_nIter = nIterStage4;
out.stage4_Pf_hat = Pf_hat_final;
out.stage4_Pf_SS = Pf_SS_final;

end

function filePath = locateRequiredFile(fileName, candidateDirs)
filePath = locateOptionalFile(fileName, candidateDirs);
assert(~isempty(filePath), 'Arquivo nao encontrado: %s (locais verificados: %s)', ...
    fileName, strjoin(candidateDirs, ', '));
end

function filePath = locateOptionalFile(fileName, candidateDirs)
filePath = '';
for i = 1:numel(candidateDirs)
    candidate = fullfile(candidateDirs{i}, fileName);
    if isfile(candidate)
        filePath = candidate;
        return;
    end
end
end

function outSS = extractSubsetSimulationOutput(S)
if isfield(S, 'outSS')
    outSS = S.outSS;
    return;
end

fields = fieldnames(S);
for i = 1:numel(fields)
    if isstruct(S.(fields{i})) && any(isfield(S.(fields{i}), {'Pf','pf'}))
        outSS = S.(fields{i});
        return;
    end
end

error('Nao foi possivel localizar a estrutura de resultado da A9.');
end

function pfRef = getStage0PfRef(T)
pfRef = getFirstTableValue(T, {'Pf_ref','PfRef','pf_ref'});
assert(~isnan(pfRef), 'Nao foi possivel localizar Pf_ref em summary_stage0.csv.');
end

function value = getMethodMetric(T, methodCandidates, valueCandidates)
value = NaN;
methodCol = findColumnName(T, {'Method','method','Metodo','Label','Name'});
valueCol = findColumnName(T, valueCandidates);

if ~isempty(methodCol) && ~isempty(valueCol)
    methods = string(T.(methodCol));
    target = normalizeTokens(methodCandidates);
    for i = 1:numel(methods)
        if any(strcmp(normalizeOne(methods(i)), target))
            value = firstNumericFromArray(T.(valueCol)(i));
            if ~isnan(value)
                return;
            end
        end
    end
end

if isempty(valueCol)
    for i = 1:numel(valueCandidates)
        valueCol = findColumnName(T, valueCandidates(i));
        if ~isempty(valueCol)
            value = firstNumericFromArray(T.(valueCol));
            if ~isnan(value)
                return;
            end
        end
    end
end
end

function [value, selectedColumn] = extractCompareMetric(T, wideCandidates, rowMethodCandidates, longValueCandidates)
selectedColumn = '';
value = NaN;

selectedColumn = findColumnName(T, wideCandidates);
if ~isempty(selectedColumn)
    methodCol = findColumnName(T, {'Method','method','Metodo','Label','Name'});
    if ~isempty(methodCol)
        methods = string(T.(methodCol));
        targetMethods = normalizeTokens(rowMethodCandidates);
        for i = 1:numel(methods)
            methodName = normalizeOne(methods(i));
            if any(strcmp(methodName, targetMethods))
                value = firstNumericFromArray(T.(selectedColumn)(i));
                if ~isnan(value)
                    return;
                end
            end
        end
    end

    value = firstNumericFromArray(T.(selectedColumn));
    if ~isnan(value)
        return;
    end
end

methodCol = findColumnName(T, {'Method','method','Metodo','Label','Name'});
valueCol = findColumnName(T, longValueCandidates);
if ~isempty(methodCol) && ~isempty(valueCol)
    methods = string(T.(methodCol));
    targetMethods = normalizeTokens(rowMethodCandidates);
    for i = 1:numel(methods)
        methodName = normalizeOne(methods(i));
        for j = 1:numel(targetMethods)
            if strcmp(methodName, targetMethods{j})
                value = firstNumericFromArray(T.(valueCol)(i));
                selectedColumn = sprintf('%s (linha %s)', valueCol, char(methods(i)));
                if ~isnan(value)
                    return;
                end
            end
        end
    end
end
end

function value = getFirstTableValue(T, candidates)
value = NaN;
col = findColumnName(T, candidates);
if isempty(col)
    return;
end
value = firstNumericFromArray(T.(col));
if isnan(value)
    raw = T.(col);
    value = string(raw(1));
end
end

function [vec, selectedColumn] = getHistoryVector(T, candidates)
selectedColumn = findColumnName(T, candidates);
vec = [];
if isempty(selectedColumn)
    return;
end
raw = T.(selectedColumn);
numericVals = nan(size(raw, 1), 1);
for i = 1:numel(numericVals)
    numericVals(i) = firstNumericFromArray(raw(i));
end
if any(~isnan(numericVals))
    vec = numericVals;
end
end

function value = getStructMetric(S, candidates)
value = NaN;
fields = fieldnames(S);
normalizedFields = normalizeTokens(fields);
normalizedCandidates = normalizeTokens(candidates);
for i = 1:numel(candidates)
    idx = find(strcmp(normalizedFields, normalizedCandidates{i}), 1, 'first');
    if ~isempty(idx)
        value = firstNumericFromArray(S.(fields{idx}));
        return;
    end
end
end

function col = findColumnName(T, candidates)
col = '';
if isempty(T) || isempty(T.Properties.VariableNames)
    return;
end

vars = T.Properties.VariableNames;
varsNorm = normalizeTokens(vars);
candNorm = normalizeTokens(candidates);

for i = 1:numel(candidates)
    exactIdx = find(strcmp(vars, candidates{i}), 1, 'first');
    if ~isempty(exactIdx)
        col = vars{exactIdx};
        return;
    end
end

for i = 1:numel(candNorm)
    idx = find(strcmp(varsNorm, candNorm{i}), 1, 'first');
    if ~isempty(idx)
        col = vars{idx};
        return;
    end
end
end

function out = normalizeTokens(values)
if isstring(values)
    values = cellstr(values(:));
elseif ischar(values)
    values = {values};
end
out = cell(size(values));
for i = 1:numel(values)
    out{i} = normalizeOne(values{i});
end
end

function out = normalizeOne(value)
out = lower(regexprep(char(string(value)), '[^a-zA-Z0-9]', ''));
end

function value = firstNumericFromArray(raw)
value = NaN;

if isnumeric(raw) || islogical(raw)
    raw = raw(:);
    idx = find(~isnan(double(raw)), 1, 'first');
    if ~isempty(idx)
        value = double(raw(idx));
    end
    return;
end

if iscell(raw)
    for i = 1:numel(raw)
        value = firstNumericFromArray(raw{i});
        if ~isnan(value)
            return;
        end
    end
    return;
end

try
    tmp = str2double(string(raw(:)));
    idx = find(~isnan(tmp), 1, 'first');
    if ~isempty(idx)
        value = tmp(idx);
    end
catch
    value = NaN;
end
end

function x = pf2beta(pf)
if isnan(pf)
    x = NaN;
    return;
end
epsv = 1e-15;
pf = min(max(pf, epsv), 1 - epsv);
x = -norminv(pf);
end

function ratio = safeDivide(a, b)
ratio = NaN(size(a));
validMask = isfinite(a) & isfinite(b) & abs(b) >= eps;
if isscalar(b)
    ratio(validMask) = a(validMask) ./ b;
    return;
end
ratio(validMask) = a(validMask) ./ b(validMask);
end

function relErrPct = computeRelativeErrorPct(absErr, pfRef)
if abs(pfRef) < eps
    if absErr < eps
        relErrPct = 0;
    else
        relErrPct = Inf;
    end
else
    relErrPct = 100 * safeDivide(absErr, abs(pfRef));
end
end

function label = classifyRelativeBias(pf, pfRef)
if isnan(pf) || isnan(pfRef)
    label = "indisponível";
elseif abs(pf - pfRef) < max(1e-12, 1e-6 * max(abs(pfRef), 1))
    label = "coincidente";
elseif pf > pfRef
    label = "superestima";
else
    label = "subestima";
end
end

function [iterOut, valuesOut] = buildPlotSeries(iterVec, valuesVec)
iterOut = [];
valuesOut = [];

if isempty(iterVec) || isempty(valuesVec)
    return;
end

n = min(numel(iterVec), numel(valuesVec));
iterOut = iterVec(1:n);
valuesOut = valuesVec(1:n);
mask = isfinite(iterOut) & isfinite(valuesOut);
iterOut = iterOut(mask);
valuesOut = valuesOut(mask);
end

function textOut = buildMetricDetail(colPf, colBeta, colCov, colErr)
parts = strings(0, 1);
if ~isempty(colPf)
    parts(end+1) = "Pf=" + string(colPf); %#ok<AGROW>
end
if ~isempty(colBeta)
    parts(end+1) = "beta=" + string(colBeta); %#ok<AGROW>
end
if ~isempty(colCov)
    parts(end+1) = "CoV=" + string(colCov); %#ok<AGROW>
end
if ~isempty(colErr)
    parts(end+1) = "erro=" + string(colErr); %#ok<AGROW>
end
if isempty(parts)
    textOut = "sem coluna identificada";
else
    textOut = strjoin(cellstr(parts), '; ');
end
end

function plotComparisonBars(labels, values, colorRGB)
x = 1:numel(values);
bar(x, values, 'FaceColor', colorRGB);
set(gca, 'XTick', x, 'XTickLabel', cellstr(labels));
grid on;
xtickangle(25);
end

function legendStrings = buildStage4Legend(hasPfHat, hasPfSS, hasR2, hasLOO)
legendStrings = {};
if hasPfHat
    legendStrings{end+1} = 'Pf hat'; %#ok<AGROW>
end
if hasPfSS
    legendStrings{end+1} = 'Pf SS'; %#ok<AGROW>
end
if hasR2
    legendStrings{end+1} = 'R2 best'; %#ok<AGROW>
end
if hasLOO
    legendStrings{end+1} = 'LOO best'; %#ok<AGROW>
end
end

function value = fallbackLabel(label)
if isempty(label)
    value = 'nao encontrado';
else
    value = char(string(label));
end
end

function out = fmtValue(value)
if isnumeric(value)
    if isnan(value)
        out = 'NaN';
    else
        out = num2str(value);
    end
else
    out = char(string(value));
end
end

function out = formatNumber(value, fmt)
if nargin < 2
    fmt = '%.6g';
end

if isnan(value)
    out = 'NaN';
elseif isinf(value)
    if value > 0
        out = 'Inf';
    else
        out = '-Inf';
    end
else
    out = sprintf(fmt, value);
end
end

function writeInterpretation(fid, methodName, pfValue, pfRef)
if isnan(pfValue)
    fprintf(fid, '   %-32s : valor indisponivel.\n', methodName);
    return;
end

absErr = abs(pfValue - pfRef);
relErr = computeRelativeErrorPct(absErr, pfRef);
if pfValue > pfRef
    relation = 'superestima';
elseif pfValue < pfRef
    relation = 'subestima';
else
    relation = 'coincide com';
end

fprintf(fid, '   %-32s : Pf = %s, %s Pf_ref, |erro abs| = %s, erro rel = %s %%.\n', ...
    methodName, formatNumber(pfValue, '%.6g'), relation, formatNumber(absErr, '%.6g'), formatNumber(relErr, '%.3f'));
end

function safeCloseFile(fid)
if fid > 0
    fclose(fid);
end
end

function safeCloseFigure(figHandle)
if ishghandle(figHandle)
    close(figHandle);
end
end

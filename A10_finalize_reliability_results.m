function out = A10_finalize_reliability_results()
% A10_FINALIZE_RELIABILITY_RESULTS
% Consolida os artefatos das etapas anteriores em uma saida final pronta
% para uso em relatorio/dissertacao, sem expor rotulos internos do tipo A6.
%
% Entradas esperadas em out_incremental:
%   - A5_pf_comparison_summary.csv
%   - A5_pf_comparison_report.txt
%   - stage4_al_history.csv
%   - stage4_report.csv
%   - pf_compare_form_mcs.csv
%   - summary_stage0.csv
%   - A9_subset_simulation_result.mat
%
% Saidas:
%   - final_reliability_summary.csv
%   - final_reliability_report.txt
%   - final_reliability_comparison.png

    clc;

    work_dir = local_resolve_work_dir();
    out_dir = fullfile(work_dir, 'out_incremental');
    if ~isfolder(out_dir)
        error('A10_finalize_reliability_results:MissingOutputDir', ...
            'Nao encontrado diretorio: %s', out_dir);
    end

    paths = struct();
    paths.a5_summary = fullfile(out_dir, 'A5_pf_comparison_summary.csv');
    paths.a5_report = fullfile(out_dir, 'A5_pf_comparison_report.txt');
    paths.stage4_history = fullfile(out_dir, 'stage4_al_history.csv');
    paths.stage4_report = fullfile(out_dir, 'stage4_report.csv');
    paths.form_mcs = fullfile(out_dir, 'pf_compare_form_mcs.csv');
    paths.stage0 = fullfile(out_dir, 'summary_stage0.csv');
    paths.a9_mat = fullfile(out_dir, 'A9_subset_simulation_result.mat');
    paths.sobol = fullfile(out_dir, 'sobol_indices.csv');
    paths.sobol_sel = fullfile(out_dir, 'sobol_selected_cumST_99.csv');

    T_a5 = local_read_optional_table(paths.a5_summary);
    T_hist = local_read_optional_table(paths.stage4_history);
    T_stage4 = local_read_optional_table(paths.stage4_report);
    T_form = local_read_optional_table(paths.form_mcs);
    T_stage0 = local_read_optional_table(paths.stage0);
    T_sobol = local_read_optional_table(paths.sobol);
    T_sobol_sel = local_read_optional_table(paths.sobol_sel);
    txtA5 = local_read_optional_text(paths.a5_report);
    outSS = local_read_optional_outss(paths.a9_mat);

    missing = local_collect_missing(paths);

    rows = {};
    rows{end+1} = local_build_rs2_row(T_a5, T_stage0); %#ok<AGROW>
    rows{end+1} = local_build_pck_row(T_a5, T_stage4, T_hist, T_form); %#ok<AGROW>
    rows{end+1} = local_build_subset_row(T_a5, outSS, T_hist); %#ok<AGROW>
    rows{end+1} = local_build_form_row(T_form); %#ok<AGROW>
    rows{end+1} = local_build_mcs_row(T_form, T_stage0, T_a5); %#ok<AGROW>
    rows{end+1} = local_build_mcs_calibrado_row(T_form); %#ok<AGROW>

    rowMask = cellfun(@(s) isstruct(s) && isfield(s, 'Pf') && ~isnan(s.Pf), rows);
    rows = rows(rowMask);
    if isempty(rows)
        error('A10_finalize_reliability_results:NoData', ...
            'Nenhuma estimativa de Pf foi encontrada nos artefatos lidos.');
    end

    nRows = numel(rows);
    Method = strings(nRows, 1);
    Pf = nan(nRows, 1);
    beta = nan(nRows, 1);
    CoV_Pf = nan(nRows, 1);
    Source = strings(nRows, 1);

    for i = 1:nRows
        Method(i) = string(rows{i}.Method);
        Pf(i) = rows{i}.Pf;
        beta(i) = rows{i}.beta;
        CoV_Pf(i) = rows{i}.CoV_Pf;
        Source(i) = string(rows{i}.Source);
    end

    T_final = table(Method, Pf, beta, CoV_Pf, Source);
    T_final = sortrows(T_final, 'Pf', 'ascend');
    T_public = removevars(T_final, 'Source');

    final_csv = fullfile(out_dir, 'final_reliability_summary.csv');
    final_txt = fullfile(out_dir, 'final_reliability_report.txt');
    final_png = fullfile(out_dir, 'final_reliability_comparison.png');

    writetable(T_public, final_csv);
    local_write_report(final_txt, T_final, T_stage4, T_hist, txtA5, missing, outSS, T_sobol, T_sobol_sel);
    local_make_figure(final_png, T_final, T_hist);

    stage4_status = local_find_status(T_stage4);
    rs2 = local_find_method(T_final, 'RS2 empírico');
    pck = local_find_method(T_final, 'PCK');
    ss = local_find_method(T_final, 'Subset Simulation');

    fprintf('\n=== Final reliability summary ===\n');
    fprintf('Diretorio de trabalho: %s\n', work_dir);
    local_print_method_line('RS2 empirico', rs2);
    local_print_method_line('PCK', pck);
    local_print_method_line('Subset Simulation', ss);
    fprintf('Status final do aprendizado ativo: %s\n', stage4_status);
    fprintf('CSV final: %s\n', final_csv);
    fprintf('Relatorio final: %s\n', final_txt);
    fprintf('Figura final: %s\n', final_png);

    out = struct();
    out.work_dir = work_dir;
    out.out_dir = out_dir;
    out.summary = T_public;
    out.summary_with_source = T_final;
    out.paths = struct('csv', final_csv, 'report', final_txt, 'figure', final_png);
    out.stage4_status = stage4_status;
    out.missing = missing;
end

function row = local_build_rs2_row(T_a5, T_stage0)
    row = local_empty_row();
    row.Method = 'RS2 empírico';
    row.Source = '';

    [pf, src] = local_find_metric_scalar(T_a5, ...
        {'Pf_RS2_empirico', 'Pf_RS2_empirico_ref', 'Pf_RS2', 'Pf_empirical_RS2', 'Pf_ref'});
    [covv, ~] = local_find_metric_scalar(T_a5, {'CoV_RS2', 'CoV_Pf_RS2', 'CoV_ref'});

    if isnan(pf)
        [pf, src2] = local_find_metric_scalar(T_stage0, {'Pf_ref', 'Pf_RS2', 'Pf_empirico', 'Pf'});
        if ~isnan(pf)
            src = src2;
        end
    end
    if isnan(covv)
        [covv, ~] = local_find_metric_scalar(T_stage0, {'CoV_ref', 'CoV_Pf_ref', 'CoV_Pf'});
    end

    if ~isnan(pf)
        row.Pf = pf;
        row.beta = local_beta_from_pf(pf);
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_build_pck_row(T_a5, T_stage4, T_hist, T_form)
    row = local_empty_row();
    row.Method = 'PCK';
    row.Source = '';

    [pf, src] = local_find_metric_scalar(T_a5, {'Pf_PCK', 'Pf_hat', 'Pf_surrogate', 'Pf_PCK_final'});
    [covv, ~] = local_find_metric_scalar(T_a5, {'CoV_PCK', 'CoV_Pf_PCK', 'CoV_hat'});

    if isnan(pf)
        [pf, src] = local_find_metric_scalar(T_stage4, {'Pf_hat', 'Pf_PCK', 'Pf'});
    end
    if isnan(pf)
        [pf, src] = local_find_last_history_value(T_hist, {'Pf_hat'});
    end
    if isnan(covv)
        [covv, ~] = local_find_metric_scalar(T_form, {'CoV_surrogate_on_RS2', 'CoV_Pf_surrogate_on_RS2'});
    end

    if ~isnan(pf)
        row.Pf = pf;
        row.beta = local_beta_from_pf(pf);
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_build_subset_row(T_a5, outSS, T_hist)
    row = local_empty_row();
    row.Method = 'Subset Simulation';
    row.Source = '';

    pf = NaN;
    beta = NaN;
    covv = NaN;
    src = '';

    if ~isempty(outSS)
        if isfield(outSS, 'Pf'), pf = outSS.Pf; end
        if isfield(outSS, 'beta'), beta = outSS.beta; end
        if isfield(outSS, 'CoV'), covv = outSS.CoV; end
        src = 'A9_subset_simulation_result.mat';
    end

    if isnan(pf)
        [pf, src] = local_find_metric_scalar(T_a5, {'Pf_SS', 'Pf_SubSetSimulation', 'Pf_SubsetSimulation'});
    end
    if isnan(beta) && ~isnan(pf)
        beta = local_beta_from_pf(pf);
    end
    if isnan(covv)
        [covv, ~] = local_find_metric_scalar(T_a5, {'CoV_SS', 'CoV_Pf_SS'});
    end
    if isnan(pf)
        [pf, src] = local_find_last_history_value(T_hist, {'Pf_SS'});
        if ~isnan(pf) && isnan(beta)
            beta = local_beta_from_pf(pf);
        end
    end

    if ~isnan(pf)
        row.Pf = pf;
        row.beta = beta;
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_build_form_row(T_form)
    row = local_empty_row();
    row.Method = 'FORM';
    row.Source = '';

    [pf, src] = local_find_metric_scalar(T_form, {'Pf_FORM'});
    [beta, ~] = local_find_metric_scalar(T_form, {'beta_FORM'});
    [covv, ~] = local_find_metric_scalar(T_form, {'CoV_FORM', 'CoV_Pf_FORM'});

    if ~isnan(pf)
        row.Pf = pf;
        if isnan(beta), beta = local_beta_from_pf(pf); end
        row.beta = beta;
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_build_mcs_row(T_form, T_stage0, T_a5)
    row = local_empty_row();
    row.Method = 'MCS';
    row.Source = '';

    [pf, src] = local_find_metric_scalar(T_form, {'Pf_MCS', 'Pf_surrogate_on_RS2'});
    [covv, ~] = local_find_metric_scalar(T_form, {'CoV_MCS', 'CoV_surrogate_on_RS2', 'CoV_Pf_surrogate_on_RS2'});

    if isnan(pf)
        [pf, src] = local_find_metric_scalar(T_a5, {'Pf_MCS'});
    end
    if isnan(pf)
        [pf, src] = local_find_metric_scalar(T_stage0, {'Pf_MCS'});
    end

    if ~isnan(pf)
        row.Pf = pf;
        row.beta = local_beta_from_pf(pf);
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_build_mcs_calibrado_row(T_form)
    row = local_empty_row();
    row.Method = 'MCS calibrado';
    row.Source = '';

    [pf, src] = local_find_metric_scalar(T_form, {'Pf_surrogate_on_RS2_calibrated'});
    [covv, ~] = local_find_metric_scalar(T_form, {'CoV_surrogate_on_RS2_calibrated', 'CoV_Pf_surrogate_on_RS2_calibrated'});

    if ~isnan(pf)
        row.Pf = pf;
        row.beta = local_beta_from_pf(pf);
        row.CoV_Pf = covv;
        row.Source = src;
    end
end

function row = local_empty_row()
    row = struct('Method', '', 'Pf', NaN, 'beta', NaN, 'CoV_Pf', NaN, 'Source', '');
end

function T = local_read_optional_table(pathFile)
    T = table();
    if ~isfile(pathFile)
        warning('A10_finalize_reliability_results:MissingFile', ...
            'Arquivo nao encontrado: %s', pathFile);
        return;
    end
    try
        T = readtable(pathFile, 'VariableNamingRule', 'preserve');
    catch ME
        warning('A10_finalize_reliability_results:ReadTableFailed', ...
            'Falha ao ler %s (%s).', pathFile, ME.message);
        T = table();
    end
end

function txt = local_read_optional_text(pathFile)
    txt = '';
    if ~isfile(pathFile)
        warning('A10_finalize_reliability_results:MissingFile', ...
            'Arquivo nao encontrado: %s', pathFile);
        return;
    end
    fid = fopen(pathFile, 'r');
    if fid < 0
        warning('A10_finalize_reliability_results:ReadTextFailed', ...
            'Nao foi possivel abrir %s.', pathFile);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    txt = fread(fid, '*char')';
    txt = strtrim(txt);
    clear cleaner;
end

function outSS = local_read_optional_outss(pathFile)
    outSS = [];
    if ~isfile(pathFile)
        warning('A10_finalize_reliability_results:MissingFile', ...
            'Arquivo nao encontrado: %s', pathFile);
        return;
    end
    S = load(pathFile);
    if isfield(S, 'outSS')
        outSS = S.outSS;
    else
        warning('A10_finalize_reliability_results:MissingOutSS', ...
            'Arquivo MAT sem struct outSS: %s', pathFile);
    end
end

function missing = local_collect_missing(paths)
    names = fieldnames(paths);
    missing = strings(0, 1);
    for i = 1:numel(names)
        if ~isfile(paths.(names{i}))
            missing(end+1, 1) = string(paths.(names{i})); %#ok<AGROW>
        end
    end
end

function [value, source] = local_find_metric_scalar(T, candidateNames)
    value = NaN;
    source = '';
    if isempty(T) || width(T) == 0 || height(T) == 0
        return;
    end

    for i = 1:numel(candidateNames)
        idx = local_find_var(T.Properties.VariableNames, candidateNames{i});
        if ~isempty(idx)
            value = local_first_numeric(T{:, idx});
            if ~isnan(value)
                source = char(T.Properties.VariableNames{idx});
                return;
            end
        end
    end

    % fallback para formato "long" com coluna de nome/metodo + valor
    textCols = local_text_columns(T);
    valueCols = local_value_columns(T, candidateNames);
    if isempty(textCols) || isempty(valueCols)
        return;
    end

    for ic = 1:numel(textCols)
        labels = local_column_strings(T{:, textCols(ic)});
        for i = 1:numel(candidateNames)
            target = local_norm_name(candidateNames{i});
            hit = contains(local_norm_name(labels), target);
            if any(hit)
                for vc = 1:numel(valueCols)
                    value = local_first_numeric(T{hit, valueCols(vc)});
                    if ~isnan(value)
                        source = sprintf('%s/%s', T.Properties.VariableNames{textCols(ic)}, T.Properties.VariableNames{valueCols(vc)});
                        return;
                    end
                end
            end
        end
    end
end

function [value, source] = local_find_last_history_value(T, candidateNames)
    value = NaN;
    source = '';
    if isempty(T) || width(T) == 0 || height(T) == 0
        return;
    end
    for i = 1:numel(candidateNames)
        idx = local_find_var(T.Properties.VariableNames, candidateNames{i});
        if ~isempty(idx)
            value = local_last_numeric(T{:, idx});
            if ~isnan(value)
                source = char(T.Properties.VariableNames{idx});
                return;
            end
        end
    end
end

function idx = local_find_var(varNames, candidate)
    idx = [];
    normVars = local_norm_name(string(varNames));
    normCandidate = local_norm_name(string(candidate));

    hit = find(normVars == normCandidate, 1, 'first');
    if ~isempty(hit)
        idx = hit;
        return;
    end

    hit = find(contains(normVars, normCandidate), 1, 'first');
    if ~isempty(hit)
        idx = hit;
    end
end

function beta = local_beta_from_pf(Pf)
    Pf = max(min(Pf, 1 - 1e-15), 1e-15);
    beta = -sqrt(2) * erfcinv(2 * Pf);
end

function s = local_norm_name(s)
    s = lower(string(s));
    s = regexprep(s, '[^a-z0-9]+', '');
end

function v = local_first_numeric(data)
    vec = local_numeric_vector(data);
    idx = find(isfinite(vec), 1, 'first');
    if isempty(idx)
        v = NaN;
    else
        v = vec(idx);
    end
end

function v = local_last_numeric(data)
    vec = local_numeric_vector(data);
    idx = find(isfinite(vec), 1, 'last');
    if isempty(idx)
        v = NaN;
    else
        v = vec(idx);
    end
end

function vec = local_numeric_vector(data)
    if isnumeric(data) || islogical(data)
        vec = double(data(:));
        return;
    end
    if iscell(data)
        vec = nan(numel(data), 1);
        for i = 1:numel(data)
            vec(i) = local_scalar_to_double(data{i});
        end
        return;
    end
    if isstring(data) || ischar(data) || iscategorical(data)
        strs = local_column_strings(data);
        vec = nan(numel(strs), 1);
        for i = 1:numel(strs)
            vec(i) = str2double(strs(i));
        end
        return;
    end
    try
        vec = double(data(:));
    catch
        vec = nan(numel(data), 1);
    end
end

function value = local_scalar_to_double(x)
    if isnumeric(x) || islogical(x)
        value = double(x(1));
    elseif isstring(x) || ischar(x)
        value = str2double(string(x));
    else
        value = NaN;
    end
end

function cols = local_text_columns(T)
    cols = [];
    for i = 1:width(T)
        col = T{:, i};
        if iscellstr(col) || isstring(col) || ischar(col) || iscategorical(col) || iscell(col)
            cols(end+1) = i; %#ok<AGROW>
        end
    end
end

function cols = local_value_columns(T, candidateNames)
    cols = [];
    for i = 1:width(T)
        nameOk = false;
        for j = 1:numel(candidateNames)
            if ~isempty(local_find_var(T.Properties.VariableNames(i), candidateNames{j}))
                nameOk = true;
                break;
            end
        end
        if nameOk || isnumeric(T{:, i}) || islogical(T{:, i})
            cols(end+1) = i; %#ok<AGROW>
        end
    end
end

function strs = local_column_strings(data)
    if iscell(data)
        strs = strings(numel(data), 1);
        for i = 1:numel(data)
            strs(i) = string(data{i});
        end
        return;
    end
    if ischar(data)
        strs = string(cellstr(data));
        return;
    end
    strs = string(data(:));
end

function local_write_report(pathFile, T_final, T_stage4, T_hist, txtA5, missing, outSS, T_sobol, T_sobol_sel)
    fid = fopen(pathFile, 'w');
    if fid < 0
        error('A10_finalize_reliability_results:WriteFailed', ...
            'Nao foi possivel gravar %s.', pathFile);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, 'FINAL RELIABILITY REPORT\n');
    fprintf(fid, 'Gerado em: %s\n\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    fprintf(fid, 'Tabela final de comparacao:\n');
    for i = 1:height(T_final)
        fprintf(fid, ' - %-20s Pf=%-12.6g beta=%-10.4f CoV(Pf)=%s\n', ...
            char(T_final.Method(i)), T_final.Pf(i), T_final.beta(i), ...
            local_num_or_na(T_final.CoV_Pf(i)));
    end
    fprintf(fid, '\n');

    stage4_status = local_find_status(T_stage4);
    fprintf(fid, 'Status final do aprendizado ativo: %s\n', stage4_status);

    [pfHat, srcPfHat] = local_find_last_history_value(T_hist, {'Pf_hat'});
    [pfSS, srcPfSS] = local_find_last_history_value(T_hist, {'Pf_SS'});
    [r2best, srcR2] = local_find_last_history_value(T_hist, {'R2_best'});
    [loobest, srcLoo] = local_find_last_history_value(T_hist, {'LOO_best'});

    fprintf(fid, 'Ultimos indicadores do historico AL:\n');
    fprintf(fid, ' - Pf_hat: %s (%s)\n', local_num_or_na(pfHat), srcPfHat);
    fprintf(fid, ' - Pf_SS : %s (%s)\n', local_num_or_na(pfSS), srcPfSS);
    fprintf(fid, ' - R2_best: %s (%s)\n', local_num_or_na(r2best), srcR2);
    fprintf(fid, ' - LOO_best: %s (%s)\n', local_num_or_na(loobest), srcLoo);

    if ~isempty(outSS)
        fprintf(fid, '\nResumo Subset Simulation (A9):\n');
        ssFields = {'Pf', 'beta', 'CoV', 'nLevels', 'pLast', 'failReached'};
        for i = 1:numel(ssFields)
            if isfield(outSS, ssFields{i})
                value = outSS.(ssFields{i});
                if isnumeric(value) || islogical(value)
                    fprintf(fid, ' - %s: %s\n', ssFields{i}, mat2str(value));
                end
            end
        end
    end

    if ~isempty(T_sobol)
        fprintf(fid, '\nResultados de sensibilidade global disponiveis.\n');
    end
    if ~isempty(T_sobol_sel)
        fprintf(fid, 'Selecao acumulada por ST>=0.99 disponivel.\n');
    end

    if strlength(string(txtA5)) > 0
        fprintf(fid, '\nResumo textual previo considerado na consolidacao final.\n');
    end

    if ~isempty(missing)
        fprintf(fid, '\nObservacao: parte dos insumos esperados nao estava disponivel, e a consolidacao utilizou apenas os resultados encontrados.\n');
        for i = 1:numel(missing)
            fprintf(fid, ' - Insumo opcional ausente %d.\n', i);
        end
    end

    clear cleaner;
end

function local_make_figure(pathFile, T_final, T_hist)
    fig = figure('Color', 'w', 'Position', [80 80 1400 900]);
    figCleaner = onCleanup(@() local_close_figure(fig));
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); %#ok<NASGU>

    nexttile;
    local_plot_method_metric(T_final.Method, T_final.Pf, 'Pf', true);
    title('Comparacao final de Pf');

    nexttile;
    local_plot_method_metric(T_final.Method, T_final.beta, '\beta', false);
    title('Comparacao final de \beta');

    nexttile;
    local_plot_pf_history(T_hist);
    title('Evolucao de Pf ao longo do AL');

    nexttile;
    local_plot_quality_history(T_hist);
    title('Qualidade do surrogate ao longo do AL');

    exportgraphics(fig, pathFile, 'Resolution', 250);
    clear figCleaner;
    local_close_figure(fig);
end

function local_plot_method_metric(methods, values, xlab, useLog)
    y = 1:numel(methods);
    good = isfinite(values);
    if useLog
        good = good & (values > 0);
    end
    if any(good)
        yPlot = y(good);
        plot(values(good), yPlot, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
            'Color', [0.1 0.35 0.75], 'MarkerFaceColor', [0.1 0.35 0.75]);
        grid on;
        set(gca, 'YDir', 'reverse');
        set(gca, 'YTick', yPlot, 'YTickLabel', string(methods(good)));
        xlabel(xlab);
        if useLog
            set(gca, 'XScale', 'log');
        end
    else
        axis off;
        text(0.5, 0.5, 'Dados indisponiveis', 'HorizontalAlignment', 'center');
    end
end

function local_plot_pf_history(T_hist)
    if isempty(T_hist) || height(T_hist) == 0
        axis off;
        text(0.5, 0.5, 'stage4_al_history.csv indisponivel', ...
            'HorizontalAlignment', 'center');
        return;
    end

    [iter, ~] = local_find_history_axis(T_hist);
    [pfHat, ~] = local_find_history_series(T_hist, {'Pf_hat'});
    [pfSS, ~] = local_find_history_series(T_hist, {'Pf_SS'});

    hold on;
    plotted = false;
    if any(isfinite(pfHat))
        semilogy(iter, pfHat, '-o', 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Pf_{hat}');
        plotted = true;
    end
    if any(isfinite(pfSS))
        semilogy(iter, pfSS, '-s', 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Pf_{SS}');
        plotted = true;
    end
    grid on;
    xlabel('Iteracao AL');
    ylabel('Pf');
    if plotted
        legend('Location', 'best');
    else
        text(0.5, 0.5, 'Series Pf_hat/Pf_SS indisponiveis', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
    end
end

function local_plot_quality_history(T_hist)
    if isempty(T_hist) || height(T_hist) == 0
        axis off;
        text(0.5, 0.5, 'stage4_al_history.csv indisponivel', ...
            'HorizontalAlignment', 'center');
        return;
    end

    [iter, ~] = local_find_history_axis(T_hist);
    [r2best, ~] = local_find_history_series(T_hist, {'R2_best'});
    [loobest, ~] = local_find_history_series(T_hist, {'LOO_best'});

    hasR2 = any(isfinite(r2best));
    hasLoo = any(isfinite(loobest));

    if hasR2
        yyaxis left;
        plot(iter, r2best, '-o', 'LineWidth', 1.4, 'MarkerSize', 5, 'Color', [0 0.45 0.74]);
        ylabel('R2_{best}');
    end
    if hasLoo
        yyaxis right;
        plot(iter, loobest, '-s', 'LineWidth', 1.4, 'MarkerSize', 5, 'Color', [0.85 0.33 0.1]);
        ylabel('LOO_{best}');
    end

    grid on;
    xlabel('Iteracao AL');

    if ~(hasR2 || hasLoo)
        text(0.5, 0.5, 'Series R2_best/LOO_best indisponiveis', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
    end
end

function [iter, src] = local_find_history_axis(T_hist)
    [iter, src] = local_find_history_series(T_hist, {'iter', 'iteration', 'AL_iter', 'k'});
    if all(~isfinite(iter))
        iter = (1:height(T_hist))';
        src = 'row_index';
    end
end

function [series, src] = local_find_history_series(T_hist, candidates)
    series = nan(height(T_hist), 1);
    src = '';
    if isempty(T_hist) || height(T_hist) == 0
        return;
    end

    for i = 1:numel(candidates)
        idx = local_find_var(T_hist.Properties.VariableNames, candidates{i});
        if ~isempty(idx)
            series = local_numeric_vector(T_hist{:, idx});
            src = char(T_hist.Properties.VariableNames{idx});
            return;
        end
    end
end

function status = local_find_status(T_stage4)
    status = 'indisponivel';
    if isempty(T_stage4) || width(T_stage4) == 0 || height(T_stage4) == 0
        return;
    end
    idx = local_find_var(T_stage4.Properties.VariableNames, 'status');
    if isempty(idx)
        idx = local_find_var(T_stage4.Properties.VariableNames, 'final_status');
    end
    if ~isempty(idx)
        vals = local_column_strings(T_stage4{:, idx});
        vals = vals(strlength(vals) > 0);
        if ~isempty(vals)
            status = char(vals(end));
        end
    end
end

function row = local_find_method(T_final, methodName)
    row = [];
    if isempty(T_final)
        return;
    end
    hit = strcmp(string(T_final.Method), string(methodName));
    if any(hit)
        row = T_final(find(hit, 1, 'first'), :);
    end
end

function local_print_method_line(label, row)
    if isempty(row)
        fprintf('%s: n/d\n', label);
        return;
    end
    fprintf('%s: Pf=%s | beta=%s | CoV=%s\n', label, ...
        local_num_or_na(row.Pf), local_num_or_na(row.beta), local_num_or_na(row.CoV_Pf));
end

function txt = local_num_or_na(x)
    if isempty(x) || ~all(isfinite(x))
        txt = 'n/d';
    else
        txt = sprintf('%.6g', x(1));
    end
end

function work_dir = local_resolve_work_dir()
    here = pwd;
    script_dir = fileparts(mfilename('fullpath'));
    candidates = {here, script_dir, 'C:\Kazuo-Script'};

    for i = 1:numel(candidates)
        base = candidates{i};
        if isempty(base) || ~isfolder(base)
            continue;
        end
        outDir = fullfile(base, 'out_incremental');
        if isfolder(outDir)
            work_dir = base;
            return;
        end
    end

    % fallback: usa a pasta do script para facilitar execucao a partir do repo.
    if isfolder(script_dir)
        work_dir = script_dir;
        return;
    end

    error('A10_finalize_reliability_results:WorkDirNotFound', ...
        'Nao foi possivel localizar o diretorio do workflow.');
end

function local_close_figure(fig)
    if ~isempty(fig) && isgraphics(fig)
        close(fig);
    end
end

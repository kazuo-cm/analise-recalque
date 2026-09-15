function out = A9_plot_subset_simulation_diagnostics(opts)
% A9_PLOT_SUBSET_SIMULATION_DIAGNOSTICS
% Gera figura diagnóstica do Subset Simulation por nível.
%
% Interpretação da figura para dissertação:
%  - Cada cor representa um nível SS-k (subset condicional).
%  - A concentração progressiva das amostras indica aproximação da região
%    de falha.
%  - A curva vermelha g=0 (quando disponível) aproxima a fronteira de falha
%    no plano 2D adotado.
%
% Requer A9_subset_simulation_result.mat com outSS.levelDiagnostics.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end

    work_dir = local_resolve_work_dir();
    out_dir = fullfile(work_dir, 'out_incremental');
    in_mat = fullfile(out_dir, 'A9_subset_simulation_result.mat');

    if ~isfile(in_mat)
        error('A9_plot_subset_simulation_diagnostics:MissingResult', ...
            'Nao encontrado: %s', in_mat);
    end

    if ~isfield(opts, 'projection') || isempty(opts.projection), opts.projection = 'auto'; end
    if ~isfield(opts, 'bgColor') || isempty(opts.bgColor), opts.bgColor = 'w'; end
    if ~isfield(opts, 'fontSize') || isempty(opts.fontSize), opts.fontSize = 12; end
    if ~isfield(opts, 'titleFontSize') || isempty(opts.titleFontSize), opts.titleFontSize = 14; end
    if ~isfield(opts, 'lineWidth') || isempty(opts.lineWidth), opts.lineWidth = 1.8; end
    if ~isfield(opts, 'markerSize') || isempty(opts.markerSize), opts.markerSize = 20; end
    if ~isfield(opts, 'dpi') || isempty(opts.dpi), opts.dpi = 260; end

    S = load(in_mat, 'outSS');
    if ~isfield(S, 'outSS')
        error('A9_plot_subset_simulation_diagnostics:InvalidResult', ...
            'Campo outSS nao encontrado em %s', in_mat);
    end
    outSS = S.outSS;
    if ~isfield(outSS, 'levelDiagnostics') || isempty(outSS.levelDiagnostics)
        error('A9_plot_subset_simulation_diagnostics:MissingDiagnostics', ...
            ['Diagnosticos por nivel nao encontrados. Rode A9 com ' ...
             'opts.recordDiagnostics=true.']);
    end

    [ZLevels, gLevels, labels, projInfo] = local_build_projection(outSS, opts);
    nLevels = numel(ZLevels);
    cmap = lines(max(3, nLevels));
    markers = {'o','s','^','d','v','>','<','p','h','x'};

    f = figure('Color', opts.bgColor, 'Position', [120 120 1180 900], 'Visible', 'off');
    ax = axes('Parent', f);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', opts.fontSize, 'LineWidth', 1.0);

    for k = 1:nLevels
        Zk = ZLevels{k};
        gk = gLevels{k};
        if isempty(Zk)
            continue;
        end
        isFail = gk <= 0;
        mk = markers{mod(k - 1, numel(markers)) + 1};
        scatter(ax, Zk(~isFail, 1), Zk(~isFail, 2), opts.markerSize, ...
            'Marker', mk, 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.30, 'MarkerEdgeAlpha', 0.55, 'DisplayName', labels{k});
        if any(isFail)
            scatter(ax, Zk(isFail, 1), Zk(isFail, 2), opts.markerSize + 12, ...
                'Marker', mk, 'MarkerFaceColor', cmap(k, :), 'MarkerEdgeColor', [0 0 0], ...
                'LineWidth', 0.8, 'DisplayName', sprintf('%s (g\\leq0)', labels{k}));
        end
    end

    if strcmpi(projInfo.mode, 'xg')
        yline(ax, 0, 'r-', 'LineWidth', opts.lineWidth, ...
            'DisplayName', 'Fronteira de falha exata (g = 0)');
        hasContour = true;
        contourMsg = '';
    else
        [hasContour, contourMsg] = local_plot_failure_contour(ax, ZLevels, gLevels, opts);
    end
    if hasContour
        contourLine = findobj(ax, 'Tag', 'failureContour');
        if ~isempty(contourLine)
            contourLine.DisplayName = 'Fronteira de falha aproximada (g = 0)';
        end
    else
        text(ax, 0.02, 0.02, contourMsg, 'Units', 'normalized', ...
            'Color', [0.70 0.10 0.10], 'FontSize', opts.fontSize - 1, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end

    xlabel(ax, projInfo.xLabel, 'FontSize', opts.fontSize + 1);
    ylabel(ax, projInfo.yLabel, 'FontSize', opts.fontSize + 1);
    title(ax, { ...
        'Subset Simulation - diagnóstico por níveis', ...
        projInfo.description}, ...
        'FontSize', opts.titleFontSize, 'FontWeight', 'bold');
    legend(ax, 'Location', 'bestoutside');

    txt = sprintf(['Leitura: SS-1, SS-2, ... mostram os subconjuntos condicionais.\n' ...
                   'A migração das amostras para g\\leq0 indica como o método alcança a falha.']);
    annotation(f, 'textbox', [0.07 0.01 0.9 0.08], 'String', txt, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontSize', opts.fontSize - 1, 'Interpreter', 'tex');

    out_png = fullfile(out_dir, 'A9_subset_simulation_diagnostics.png');
    exportgraphics(f, out_png, 'Resolution', opts.dpi);
    close(f);

    out = struct();
    out.file_png = out_png;
    out.projection = projInfo;
    out.hasFailureContour = hasContour;
    out.contourMessage = contourMsg;
    fprintf('\n=== A9 subset simulation diagnostics ===\n');
    fprintf('Figura: %s\n', out_png);
    fprintf('Projecao: %s\n', projInfo.description);
    fprintf('Fronteira g=0: %d\n', hasContour);
end

function [ZLevels, gLevels, labels, projInfo] = local_build_projection(outSS, opts)
    LD = outSS.levelDiagnostics;
    nLevels = numel(LD);
    M = size(LD(1).X, 2);

    Xcells = cell(nLevels, 1);
    gcells = cell(nLevels, 1);
    for k = 1:nLevels
        Xcells{k} = LD(k).X;
        gcells{k} = LD(k).g(:);
    end
    Xall = vertcat(Xcells{:});
    gall = vertcat(gcells{:});

    if M == 1
        ZLevels = cell(nLevels, 1);
        gLevels = cell(nLevels, 1);
        labels = cell(nLevels, 1);
        for k = 1:nLevels
            xk = LD(k).X(:, 1);
            ZLevels{k} = [xk, LD(k).g(:)];
            gLevels{k} = LD(k).g(:);
            labels{k} = sprintf('SS-%d (b=%.3g)', LD(k).level, LD(k).threshold);
        end
        projInfo = struct();
        projInfo.mode = 'xg';
        projInfo.usesResponseOnY = true;
        projInfo.variables = 1;
        projInfo.xLabel = 'X_1';
        projInfo.yLabel = 'g(X)';
        projInfo.description = 'Modelo 1D: plano [X_1, g(X)]';
        projInfo.coordinateMeaning = 'col1 = X_1, col2 = g(X)';
        return;
    end

    modeReq = lower(string(opts.projection));
    supportedModes = ["auto","pca","influential"];
    if ~any(modeReq == supportedModes)
        error('A9_plot_subset_simulation_diagnostics:InvalidProjection', ...
            'opts.projection="%s" invalido. Use: auto, influential ou pca.', char(modeReq));
    end
    usePCA = false;
    if modeReq == "pca"
        usePCA = true;
    elseif modeReq == "auto"
        usePCA = (M > 2);
    end

    if ~usePCA && M >= 2
        if M == 2
            P = [1 2];
        else
            idx = local_pick_most_influential_dims(Xall, gall);
            P = idx(1:2);
        end
        [ZLevels, gLevels, labels, projInfo] = local_make_influential_projection(LD, P, 'influential');
        return;
    end

    finiteRows = all(isfinite(Xall), 2);
    Xpca = Xall(finiteRows, :);
    hasPCA = (exist('pca', 'file') == 2) || (exist('pca', 'builtin') == 5);
    if ~hasPCA
        idx = local_pick_most_influential_dims(Xall, gall);
        [ZLevels, gLevels, labels, projInfo] = local_make_influential_projection(LD, idx(1:2), 'influential-fallback');
        projInfo.description = sprintf(['PCA indisponível, projeção de fallback ' ...
            'nas variáveis X_%d e X_%d'], projInfo.variables(1), projInfo.variables(2));
        return;
    end
    if size(Xpca, 1) < 2
        idx = local_pick_most_influential_dims(Xall, gall);
        [ZLevels, gLevels, labels, projInfo] = local_make_influential_projection(LD, idx(1:2), 'influential-fallback');
        projInfo.description = sprintf(['PCA sem amostras finitas suficientes, projeção de fallback ' ...
            'nas variáveis X_%d e X_%d'], projInfo.variables(1), projInfo.variables(2));
        return;
    end

    [Xstd, mu, sg] = local_standardize(Xpca);
    [coeff, score] = pca(Xstd, 'NumComponents', 2);
    if size(coeff, 2) < 2 || size(score, 2) < 2
        idx = local_pick_most_influential_dims(Xall, gall);
        [ZLevels, gLevels, labels, projInfo] = local_make_influential_projection(LD, idx(1:2), 'influential-fallback');
        projInfo.description = sprintf(['PCA com posto insuficiente, projeção de fallback ' ...
            'nas variáveis X_%d e X_%d'], projInfo.variables(1), projInfo.variables(2));
        return;
    end
    ZLevels = cell(nLevels, 1);
    gLevels = cell(nLevels, 1);
    labels = cell(nLevels, 1);
    for k = 1:nLevels
        Xk = LD(k).X;
        XkStd = bsxfun(@rdivide, bsxfun(@minus, Xk, mu), sg);
        bad = any(~isfinite(XkStd), 2);
        XkStd(bad, :) = NaN;
        ZLevels{k} = XkStd * coeff(:, 1:2);
        gLevels{k} = LD(k).g(:);
        labels{k} = sprintf('SS-%d (b=%.3g)', LD(k).level, LD(k).threshold);
    end
    expVar = 100 * var(score(:, 1:2), 0, 1) ./ max(sum(var(Xstd, 0, 1)), eps);
    projInfo = struct();
    projInfo.mode = 'pca';
    projInfo.usesResponseOnY = false;
    projInfo.variables = [1 2];
    projInfo.xLabel = sprintf('PC1 (%.1f%%)', expVar(1));
    projInfo.yLabel = sprintf('PC2 (%.1f%%)', expVar(2));
    projInfo.description = sprintf('Projeção PCA 2D (componentes principais de X), var.=%.1f%%', sum(expVar));
    projInfo.coordinateMeaning = 'col1 e col2 = coordenadas PCA de X';
    projInfo.loadings = coeff;
end

function [Xs, mu, sg] = local_standardize(X)
    mu = mean(X, 1, 'omitnan');
    sg = std(X, 0, 1, 'omitnan');
    sg(~isfinite(sg) | sg <= 0) = 1;
    Xs = bsxfun(@rdivide, bsxfun(@minus, X, mu), sg);
    Xs(~isfinite(Xs)) = 0;
end

function [ZLevels, gLevels, labels, projInfo] = local_make_influential_projection(LD, P, modeName)
    nLevels = numel(LD);
    ZLevels = cell(nLevels, 1);
    gLevels = cell(nLevels, 1);
    labels = cell(nLevels, 1);
    for k = 1:nLevels
        ZLevels{k} = LD(k).X(:, P);
        gLevels{k} = LD(k).g(:);
        labels{k} = sprintf('SS-%d (b=%.3g)', LD(k).level, LD(k).threshold);
    end
    projInfo = struct();
    projInfo.mode = char(modeName);
    projInfo.usesResponseOnY = false;
    projInfo.variables = P;
    projInfo.xLabel = sprintf('X_{%d}', P(1));
    projInfo.yLabel = sprintf('X_{%d}', P(2));
    projInfo.description = sprintf('Projeção nas variáveis mais influentes: X_%d e X_%d', P(1), P(2));
    projInfo.coordinateMeaning = 'col1 e col2 = coordenadas projetadas de X';
end

function idx = local_pick_most_influential_dims(X, g)
    M = size(X, 2);
    s = zeros(1, M);
    g = g(:);
    for j = 1:M
        xj = X(:, j);
        ok = isfinite(xj) & isfinite(g);
        if nnz(ok) < 3
            s(j) = 0;
            continue;
        end
        C = corrcoef(xj(ok), g(ok));
        if numel(C) < 4
            cj = 0;
        else
            cj = C(1, 2);
            if ~isfinite(cj)
                cj = 0;
            end
        end
        s(j) = abs(cj);
    end
    [~, ord] = sort(s, 'descend');
    if M == 1
        idx = [1 1];
    else
        idx = ord(1:2);
    end
end

function [ok, msg] = local_plot_failure_contour(ax, ZLevels, gLevels, opts)
    Zcells = cell(numel(ZLevels), 1);
    gcells = cell(numel(gLevels), 1);
    for k = 1:numel(ZLevels)
        Zcells{k} = ZLevels{k};
        gcells{k} = gLevels{k}(:);
    end
    Z = vertcat(Zcells{:});
    g = vertcat(gcells{:});
    ok = false;
    msg = 'Fronteira g=0 indisponível: exibindo somente dispersão por níveis.';

    is2d = (size(Z, 2) == 2);
    hasBothSigns = any(g <= 0) && any(g > 0);
    if ~is2d || ~hasBothSigns
        return;
    end

    okFinite = isfinite(g) & isfinite(Z(:, 1)) & isfinite(Z(:, 2));
    Z = Z(okFinite, :);
    g = g(okFinite);
    Zr = round(Z * 1e8) / 1e8;
    if size(unique(Zr, 'rows'), 1) < 10
        return;
    end

    try
        F = scatteredInterpolant(Z(:, 1), Z(:, 2), g, 'natural', 'none');
        pad = 0.05;
        xlimData = [min(Z(:, 1)), max(Z(:, 1))];
        ylimData = [min(Z(:, 2)), max(Z(:, 2))];
        dx = max(diff(xlimData), eps);
        dy = max(diff(ylimData), eps);
        xx = linspace(xlimData(1) - pad * dx, xlimData(2) + pad * dx, 120);
        yy = linspace(ylimData(1) - pad * dy, ylimData(2) + pad * dy, 120);
        [XX, YY] = meshgrid(xx, yy);
        GG = F(XX, YY);
        if all(~isfinite(GG(:)))
            return;
        end
        [~, h] = contour(ax, XX, YY, GG, [0 0], 'r-', 'LineWidth', opts.lineWidth);
        set(h, 'Tag', 'failureContour');
        ok = true;
        msg = '';
    catch
        ok = false;
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
        if isfolder(fullfile(base, 'out_incremental'))
            work_dir = base;
            return;
        end
    end
    error('A9_plot_subset_simulation_diagnostics:WorkDirNotFound', ...
        'Nao foi possivel localizar o diretorio do workflow (out_incremental).');
end

function fig = A9_plot_subset_simulation_diagnostics(opts)
% A9_PLOT_SUBSET_SIMULATION_DIAGNOSTICS
% -------------------------------------------------------------------------
% Reads A9_subset_simulation_levels.mat and produces a publication-style
% 2D diagnostic plot for Subset Simulation.
%
% This figure is meant to explain, in dissertation language, how SS
% progressively concentrates samples toward the failure region:
%   - SS-1, SS-2, ... are the level-by-level sample clouds
%   - each b_k line/label is an intermediate threshold g(x) <= b_k
%   - the contour g(x) = 0 is the failure boundary
%
% Usage from the project root:
%   A9_plot_subset_simulation_diagnostics();
%   A9_plot_subset_simulation_diagnostics(struct( ...
%       'backgroundColor', [1 1 1], ...
%       'fontSize', 12, ...
%       'titleFontSize', 16, ...
%       'markerSize', 28, ...
%       'lineWidth', 1.8));
%
% Optional fields in opts
%   levelsFile        : defaults to C:\Kazuo-Script\out_incremental\A9_subset_simulation_levels.mat
%   outputPng         : optional PNG export path
%   projectionMethod  : 'auto' | 'pair' | 'pca'
%   variablePair      : [i j] or {'x1','x2'} when projectionMethod='pair'
%   backgroundColor   : figure background
%   axesColor         : axes background
%   fontName          : default 'Arial'
%   fontSize          : default 12
%   titleFontSize     : default 16
%   labelFontSize     : default 13
%   legendFontSize    : default 11
%   markerSize        : default 28
%   lineWidth         : default 1.8
%   gridSize          : default 160
% -------------------------------------------------------------------------

if nargin < 1 || isempty(opts)
    opts = struct();
end
opts = applyPlotDefaults(opts);

S = load(opts.levelsFile, 'subsetLevels');
subsetLevels = S.subsetLevels;
[projectedLevels, axisNames, basisLabel, basisMethod] = resolveProjectedLevels(subsetLevels, opts);

fig = figure('Color', opts.backgroundColor, 'Name', 'Subset Simulation diagnostics');
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'Color', opts.axesColor, ...
    'FontName', opts.fontName, ...
    'FontSize', opts.fontSize, ...
    'LineWidth', 1.0);

colors = lines(max(subsetLevels.levelCount, 1));
legendHandles = gobjects(subsetLevels.levelCount + 1, 1);
legendLabels = strings(subsetLevels.levelCount + 1, 1);

for k = 1:subsetLevels.levelCount
    Z = projectedLevels(k).samples2D;
    legendHandles(k) = scatter(ax, Z(:, 1), Z(:, 2), opts.markerSize, ...
        'MarkerFaceColor', colors(k, :), ...
        'MarkerEdgeColor', colors(k, :), ...
        'MarkerFaceAlpha', 0.30, ...
        'MarkerEdgeAlpha', 0.50, ...
        'DisplayName', sprintf('SS-%d', k));

    legendLabels(k) = sprintf('SS-%d  (b_%d = %.4g)', ...
        k, k, projectedLevels(k).threshold);
end

[xGrid, yGrid, gGrid] = buildFailureContour(projectedLevels, opts.gridSize);
finiteMask = isfinite(gGrid);
if ~isempty(gGrid) && any(finiteMask(:)) ...
        && any(gGrid(finiteMask) <= 0) && any(gGrid(finiteMask) > 0)
    if strcmpi(basisMethod, 'pca')
        contourLabel = 'Projected zero-level estimate';
    else
        contourLabel = 'Failure boundary  g(x)=0';
    end
    contourHandle = contour(ax, xGrid, yGrid, gGrid, [0 0], ...
        'Color', [0.85 0.10 0.10], ...
        'LineWidth', opts.lineWidth, ...
        'DisplayName', contourLabel);
    legendHandles(subsetLevels.levelCount + 1) = contourHandle;
    legendLabels(subsetLevels.levelCount + 1) = string(contourLabel);
end

xlabel(ax, axisNames(1), 'FontName', opts.fontName, 'FontSize', opts.labelFontSize);
ylabel(ax, axisNames(2), 'FontName', opts.fontName, 'FontSize', opts.labelFontSize);
title(ax, { ...
    sprintf('Subset Simulation diagnostic plot (%s)', basisLabel), ...
    sprintf('P_f = %.4e,  \\beta = %.4f,  CoV = %.4f', subsetLevels.Pf, subsetLevels.beta, subsetLevels.CoV)}, ...
    'FontName', opts.fontName, ...
    'FontSize', opts.titleFontSize, ...
    'FontWeight', 'bold');

legendIdx = find(isgraphics(legendHandles) & strlength(legendLabels) > 0);
legend(ax, legendHandles(legendIdx), legendLabels(legendIdx), ...
    'Location', 'bestoutside', ...
    'FontName', opts.fontName, ...
    'FontSize', opts.legendFontSize);

annotation(fig, 'textbox', [0.12 0.01 0.76 0.08], ...
    'String', ['Interpretation: each cloud SS-k shows where the algorithm sampled at level k. ' ...
               'As k increases, the samples move toward the rare-event region, and the red contour ' ...
               contourExplanation(basisMethod)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontName', opts.fontName, ...
    'FontSize', opts.fontSize, ...
    'Color', [0.15 0.15 0.15]);

axis(ax, 'tight');
if ~isempty(opts.outputPng)
    outFolder = fileparts(opts.outputPng);
    if ~isempty(outFolder) && ~isfolder(outFolder)
        mkdir(outFolder);
    end
    exportgraphics(fig, opts.outputPng, 'Resolution', 220);
end
end

function opts = applyPlotDefaults(opts)
if ~isfield(opts, 'workDir') || isempty(opts.workDir)
    opts.workDir = 'C:\Kazuo-Script';
end
if ~isfield(opts, 'levelsFile') || isempty(opts.levelsFile)
    opts.levelsFile = fullfile(opts.workDir, 'out_incremental', 'A9_subset_simulation_levels.mat');
end
if ~isfield(opts, 'outputPng')
    opts.outputPng = '';
end
if ~isfield(opts, 'projectionMethod') || isempty(opts.projectionMethod)
    opts.projectionMethod = 'auto';
end
if ~isfield(opts, 'variablePair')
    opts.variablePair = [];
end
if ~isfield(opts, 'backgroundColor') || isempty(opts.backgroundColor)
    opts.backgroundColor = 'w';
end
if ~isfield(opts, 'axesColor') || isempty(opts.axesColor)
    opts.axesColor = 'w';
end
if ~isfield(opts, 'fontName') || isempty(opts.fontName)
    opts.fontName = 'Arial';
end
if ~isfield(opts, 'fontSize') || isempty(opts.fontSize)
    opts.fontSize = 12;
end
if ~isfield(opts, 'titleFontSize') || isempty(opts.titleFontSize)
    opts.titleFontSize = 16;
end
if ~isfield(opts, 'labelFontSize') || isempty(opts.labelFontSize)
    opts.labelFontSize = 13;
end
if ~isfield(opts, 'legendFontSize') || isempty(opts.legendFontSize)
    opts.legendFontSize = 11;
end
if ~isfield(opts, 'markerSize') || isempty(opts.markerSize)
    opts.markerSize = 28;
end
if ~isfield(opts, 'lineWidth') || isempty(opts.lineWidth)
    opts.lineWidth = 1.8;
end
if ~isfield(opts, 'gridSize') || isempty(opts.gridSize)
    opts.gridSize = 160;
end
end

function [levelsOut, axisNames, basisLabel, basisMethod] = resolveProjectedLevels(subsetLevels, opts)
levelsOut = subsetLevels.levels;
method = lower(string(opts.projectionMethod));

switch method
    case "auto"
        if hasStoredAutoProjection(subsetLevels)
            for i = 1:numel(levelsOut)
                levelsOut(i).samples2D = subsetLevels.levels(i).samples2D;
            end
            axisNames = string(subsetLevels.plotBasis.varNames(:));
            basisLabel = subsetLevels.plotBasis.label;
            basisMethod = char(subsetLevels.plotBasis.method);
        else
            [levelsOut, axisNames, basisLabel, basisMethod] = rebuildAutoProjection(levelsOut, subsetLevels);
        end
    case "pair"
        pair = parseVariablePair(opts.variablePair, subsetLevels.varNames);
        for i = 1:numel(levelsOut)
            Xi = subsetLevels.levels(i).samplesX;
            if pair(1) == pair(2)
                levelsOut(i).samples2D = [Xi(:, pair(1)), zeros(size(Xi, 1), 1)];
            else
                levelsOut(i).samples2D = Xi(:, pair);
            end
        end
        axisNames = string(subsetLevels.varNames(pair));
        basisLabel = sprintf('%s vs %s', axisNames(1), axisNames(2));
        basisMethod = 'pair';
    case "pca"
        [coeff, mu, sigma] = computePcaBasis(subsetLevels);
        for i = 1:numel(levelsOut)
            Xi = subsetLevels.levels(i).samplesX;
            Xs = (Xi - mu) ./ sigma;
            Zi = Xs * coeff;
            if size(Zi, 2) < 2
                Zi(:, 2) = 0;
            end
            levelsOut(i).samples2D = Zi(:, 1:2);
        end
        axisNames = ["PC1"; "PC2"];
        basisLabel = 'PCA projection (PC1 vs PC2)';
        basisMethod = 'pca';
    otherwise
        error('A9_plot_subset_simulation_diagnostics:InvalidProjection', ...
            'projectionMethod must be auto, pair or pca.');
end
end

function txt = contourExplanation(basisMethod)
if strcmpi(basisMethod, 'pca')
    txt = 'shows a projected zero-level estimate in the displayed PCA basis.';
else
    txt = 'represents the estimated failure boundary g(x)=0 in the displayed 2D basis.';
end
end

function tf = hasStoredAutoProjection(subsetLevels)
tf = isfield(subsetLevels, 'plotBasis') ...
    && isstruct(subsetLevels.plotBasis) ...
    && isfield(subsetLevels.plotBasis, 'method') ...
    && isfield(subsetLevels.plotBasis, 'varNames') ...
    && isfield(subsetLevels.plotBasis, 'label');
if ~tf
    return;
end
tf = all(arrayfun(@(lvl) isfield(lvl, 'samples2D') && size(lvl.samples2D, 2) == 2, subsetLevels.levels));
end

function [levelsOut, axisNames, basisLabel, basisMethod] = rebuildAutoProjection(levelsOut, subsetLevels)
nVars = size(levelsOut(1).samplesX, 2);
if isfield(subsetLevels, 'metadata') && isfield(subsetLevels.metadata, 'recommendedVariablePair')
    pair = subsetLevels.metadata.recommendedVariablePair;
elseif nVars == 1
    pair = [1 1];
else
    pair = [1 2];
end

for i = 1:numel(levelsOut)
    Xi = levelsOut(i).samplesX;
    if pair(1) == pair(2)
        levelsOut(i).samples2D = [Xi(:, pair(1)), zeros(size(Xi, 1), 1)];
    else
        levelsOut(i).samples2D = Xi(:, pair);
    end
end

if pair(1) == pair(2)
    axisNames = string(subsetLevels.varNames(pair));
else
    axisNames = string(subsetLevels.varNames(pair));
end
basisLabel = sprintf('%s vs %s', axisNames(1), axisNames(2));
basisMethod = 'pair';
end

function [xGrid, yGrid, gGrid] = buildFailureContour(levelsOut, gridSize)
zCells = arrayfun(@(lvl) lvl.samples2D, levelsOut, 'UniformOutput', false);
gCells = arrayfun(@(lvl) lvl.g(:), levelsOut, 'UniformOutput', false);
zAll = vertcat(zCells{:});
gAll = vertcat(gCells{:});

if size(zAll, 1) < 3 || numel(unique(zAll(:, 1))) < 2 || numel(unique(zAll(:, 2))) < 2
    xGrid = [];
    yGrid = [];
    gGrid = [];
    return;
end

xMin = min(zAll(:, 1));
xMax = max(zAll(:, 1));
yMin = min(zAll(:, 2));
yMax = max(zAll(:, 2));

xPad = 0.05 * max(xMax - xMin, 1);
yPad = 0.05 * max(yMax - yMin, 1);

xv = linspace(xMin - xPad, xMax + xPad, gridSize);
yv = linspace(yMin - yPad, yMax + yPad, gridSize);
[xGrid, yGrid] = meshgrid(xv, yv);

F = scatteredInterpolant(zAll(:, 1), zAll(:, 2), gAll, 'natural', 'none');
gGrid = F(xGrid, yGrid);
end

function pair = parseVariablePair(variablePair, varNames)
if ~exist('variablePair', 'var') || isempty(variablePair)
    error('A9_plot_subset_simulation_diagnostics:MissingVariablePair', ...
        'opts.variablePair is required when projectionMethod = ''pair''.');
end

if isnumeric(variablePair)
    pair = reshape(variablePair, 1, []);
else
    pair = zeros(1, numel(variablePair));
    for i = 1:numel(variablePair)
        matchIdx = find(string(varNames) == string(variablePair{i}), 1, 'first');
        if isempty(matchIdx)
            pair(i) = NaN;
        else
            pair(i) = matchIdx;
        end
    end
end

if numel(pair) ~= 2 || any(isnan(pair)) || any(pair < 1) || any(pair > numel(varNames))
    error('A9_plot_subset_simulation_diagnostics:InvalidVariablePair', ...
        'opts.variablePair must identify exactly two valid variables.');
end
pair = double(pair);
end

function [coeff, mu, sigma] = computePcaBasis(subsetLevels)
X = [];
for i = 1:subsetLevels.levelCount
    X = [X; subsetLevels.levels(i).samplesX]; %#ok<AGROW>
end
mu = mean(X, 1);
sigma = std(X, 0, 1);
sigma(sigma == 0) = 1;
Xs = (X - mu) ./ sigma;
[~, ~, V] = svd(Xs, 'econ');
coeff = V(:, 1:min(2, size(V, 2)));
if size(coeff, 2) < 2
    coeff(:, 2) = 0;
end
end

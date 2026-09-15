function outSS = A9_subset_simulation_pf(gfun, sampleFcn, opts)
% A9_SUBSET_SIMULATION_PF
% -------------------------------------------------------------------------
% Backward-compatible Subset Simulation (SS) driver with per-level
% diagnostic storage for publication-style figures.
%
% The reliability outputs remain the usual:
%   outSS.Pf
%   outSS.beta
%   outSS.CoV
%
% In addition, this function saves a second artifact with the level-by-level
% sample clouds, g-values, thresholds and projection metadata needed to
% reconstruct the classic SS diagnostic figure later:
%   C:\Kazuo-Script\out_incremental\A9_subset_simulation_levels.mat
%
% Minimal usage from the project root:
%   gfun      = @(X) 3.0 - X(:,1) - 0.8 .* X(:,2);
%   sampleFcn = @(n) randn(n, 2);
%   outSS = A9_subset_simulation_pf(gfun, sampleFcn);
%
% Inputs
%   gfun      : function handle, g = gfun(X), failure defined by g <= 0
%   sampleFcn : function handle, X = sampleFcn(n), independent base samples
%   opts      : optional struct
%
% Important options
%   opts.N                : samples per level (default 1000)
%   opts.p0               : conditional probability per level (default 0.1)
%   opts.maxLevels        : maximum number of SS levels (default 10)
%   opts.proposalSigma    : scalar or 1xd RW proposal scale in X-space
%   opts.varNames         : variable names for plotting
%   opts.workDir          : defaults to 'C:\Kazuo-Script'
%   opts.outDir           : defaults to fullfile(workDir,'out_incremental')
%   opts.plotBasis        : struct controlling the stored 2D basis
%                           .method = 'auto' | 'pair' | 'pca'
%                           .variablePair = [i j] or {'x1','x2'}
%
% The saved diagnostics are intentionally self-contained so A10 can ignore
% them unless a later consolidation step wants to reference the file path.
% -------------------------------------------------------------------------

if nargin < 1 || isempty(gfun) || ~isa(gfun, 'function_handle')
    error('A9_subset_simulation_pf:InvalidInput', ...
        'gfun must be a function handle with failure defined by g(x) <= 0.');
end

if nargin < 2 || isempty(sampleFcn) || ~isa(sampleFcn, 'function_handle')
    error('A9_subset_simulation_pf:InvalidInput', ...
        'sampleFcn must be a function handle of the form X = sampleFcn(n).');
end

if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = applyDefaults(opts);
validateattributes(opts.N, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(opts.maxLevels, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(opts.p0, {'numeric'}, {'scalar', '>', 0, '<', 1});
nPerLevel = opts.N;
nSeeds = max(1, round(opts.p0 * nPerLevel));
nPerChain = max(1, ceil(nPerLevel / nSeeds));

if abs(nSeeds / nPerLevel - opts.p0) > 1e-12
    warning('A9_subset_simulation_pf:RoundedP0', ...
        'opts.p0*N was rounded to %d seeds to maintain a finite number of chains.', nSeeds);
end

X = ensure2D(sampleFcn(nPerLevel));
if size(X, 1) ~= nPerLevel
    error('A9_subset_simulation_pf:InvalidSampleCount', ...
        'sampleFcn(%d) must return exactly %d samples.', nPerLevel, nPerLevel);
end
nVars = size(X, 2);
varNames = resolveVarNames(opts, nVars);
g = evaluateLimitState(gfun, X);

levelRecords = repmat(emptyLevelRecord(nVars), opts.maxLevels, 1);
thresholds = nan(opts.maxLevels, 1);
levelMeta = repmat(struct(), opts.maxLevels, 1);

for levelIdx = 1:opts.maxLevels
    [threshold, seedIdx] = computeIntermediateThreshold(g, nSeeds);
    thresholds(levelIdx) = threshold;

    levelRecords(levelIdx) = makeLevelRecord(levelIdx, X, g, threshold, seedIdx, levelMeta(levelIdx));

    if threshold <= 0
        finalLevel = levelIdx;
        break;
    end

    seedX = X(seedIdx, :);
    seedG = g(seedIdx);
    proposalSigma = resolveProposalSigma(opts, seedX);
    [X, g, meta] = conditionalLevelSamples(gfun, seedX, seedG, threshold, nPerChain, nPerLevel, proposalSigma);
    if levelIdx < opts.maxLevels
        levelMeta(levelIdx + 1) = meta;
    end
end

if ~exist('finalLevel', 'var')
    finalLevel = opts.maxLevels;
end

levelRecords = levelRecords(1:finalLevel);
thresholds = thresholds(1:finalLevel);

if isempty(levelRecords(finalLevel).samplesX)
    levelRecords(finalLevel) = makeLevelRecord(finalLevel, X, g, thresholds(finalLevel), [], []);
end

pfLast = mean(levelRecords(finalLevel).g <= 0);
Pf = (opts.p0 ^ max(finalLevel - 1, 0)) * pfLast;
beta = pf2beta(Pf);
CoV = subsetSimulationCoV(opts.p0, nPerLevel, finalLevel, pfLast);

subsetLevels = buildDiagnosticArtifact(levelRecords, thresholds, Pf, beta, CoV, opts, varNames);
subsetLevels = attachStored2DBasis(subsetLevels, opts);

outSS = struct();
outSS.Pf = Pf;
outSS.beta = beta;
outSS.CoV = CoV;
outSS.nLevels = finalLevel;
outSS.thresholds = thresholds;
outSS.levelsFile = fullfile(opts.outDir, 'A9_subset_simulation_levels.mat');
outSS.resultFile = fullfile(opts.outDir, 'A9_subset_simulation_result.mat');
outSS.diagnosticBasis = subsetLevels.plotBasis;

if ~isfolder(opts.outDir)
    mkdir(opts.outDir);
end

save(outSS.resultFile, 'outSS');
save(outSS.levelsFile, 'subsetLevels');
end

function opts = applyDefaults(opts)
if ~isfield(opts, 'N') || isempty(opts.N)
    opts.N = 1000;
end
if ~isfield(opts, 'p0') || isempty(opts.p0)
    opts.p0 = 0.1;
end
if ~isfield(opts, 'maxLevels') || isempty(opts.maxLevels)
    opts.maxLevels = 10;
end
if ~isfield(opts, 'workDir') || isempty(opts.workDir)
    opts.workDir = 'C:\Kazuo-Script';
end
if ~isfield(opts, 'outDir') || isempty(opts.outDir)
    opts.outDir = fullfile(opts.workDir, 'out_incremental');
end
if ~isfield(opts, 'plotBasis') || isempty(opts.plotBasis)
    opts.plotBasis = struct('method', 'auto');
end
end

function X = ensure2D(X)
if isvector(X)
    X = X(:);
end
if ~isnumeric(X) || isempty(X)
    error('A9_subset_simulation_pf:InvalidSamples', ...
        'sampleFcn must return a non-empty numeric matrix.');
end
end

function g = evaluateLimitState(gfun, X)
g = gfun(X);
if isrow(g)
    g = g.';
end
if ~isnumeric(g) || numel(g) ~= size(X, 1)
    error('A9_subset_simulation_pf:InvalidLimitState', ...
        'gfun must return one numeric g-value per sample.');
end
g = g(:);
end

function varNames = resolveVarNames(opts, nVars)
if isfield(opts, 'varNames') && ~isempty(opts.varNames)
    varNames = string(opts.varNames(:));
else
    varNames = "X" + string(1:nVars);
end
if numel(varNames) ~= nVars
    error('A9_subset_simulation_pf:InvalidVarNames', ...
        'opts.varNames must have one entry per variable.');
end
end

function [threshold, seedIdx] = computeIntermediateThreshold(g, nSeeds)
[gSorted, order] = sort(g, 'ascend');
seedIdx = order(1:nSeeds);
threshold = gSorted(nSeeds);
end

function proposalSigma = resolveProposalSigma(opts, seedX)
if isfield(opts, 'proposalSigma') && ~isempty(opts.proposalSigma)
    proposalSigma = opts.proposalSigma;
else
    proposalSigma = 0.20 .* max(std(seedX, 0, 1), 1e-6);
end
if isscalar(proposalSigma)
    proposalSigma = repmat(proposalSigma, 1, size(seedX, 2));
end
proposalSigma = reshape(proposalSigma, 1, []);
end

function [Xnew, gnew, meta] = conditionalLevelSamples(gfun, seedX, seedG, threshold, nPerChain, nTarget, proposalSigma)
nSeeds = size(seedX, 1);
nVars = size(seedX, 2);
nAlloc = nSeeds * nPerChain;

Xnew = zeros(nAlloc, nVars);
gnew = zeros(nAlloc, 1);
chainId = zeros(nAlloc, 1);
accepted = false(nAlloc, 1);
isSeed = false(nAlloc, 1);

cursor = 0;
for i = 1:nSeeds
    xCurr = seedX(i, :);
    gCurr = seedG(i);
    for j = 1:nPerChain
        cursor = cursor + 1;
        if j == 1
            isSeed(cursor) = true;
            accepted(cursor) = true;
        else
            xProp = xCurr + proposalSigma .* randn(1, nVars);
            gProp = evaluateLimitState(gfun, xProp);
            if gProp <= threshold
                xCurr = xProp;
                gCurr = gProp;
                accepted(cursor) = true;
            end
        end
        Xnew(cursor, :) = xCurr;
        gnew(cursor) = gCurr;
        chainId(cursor) = i;
    end
end

Xnew = Xnew(1:nTarget, :);
gnew = gnew(1:nTarget);
chainId = chainId(1:nTarget);
accepted = accepted(1:nTarget);
isSeed = isSeed(1:nTarget);

meta = struct();
meta.chainId = chainId;
meta.accepted = accepted;
meta.isSeed = isSeed;
meta.thresholdUsed = threshold;
meta.proposalSigma = proposalSigma;
end

function level = emptyLevelRecord(nVars)
level = struct( ...
    'levelIndex', [], ...
    'samplesX', zeros(0, nVars), ...
    'samples', zeros(0, nVars), ...
    'g', zeros(0, 1), ...
    'gValues', zeros(0, 1), ...
    'threshold', NaN, ...
    'seedIndex', zeros(0, 1), ...
    'generationMeta', struct(), ...
    'samples2D', zeros(0, 2), ...
    'projectionLabel', "" ...
    );
end

function level = makeLevelRecord(levelIndex, X, g, threshold, seedIdx, generationMeta)
if nargin < 5 || isempty(seedIdx)
    seedIdx = zeros(0, 1);
end
if nargin < 6 || isempty(generationMeta)
    generationMeta = struct();
end

level = struct();
level.levelIndex = levelIndex;
level.samplesX = X;
level.samples = X;
level.g = g(:);
level.gValues = g(:);
level.threshold = threshold;
level.seedIndex = seedIdx(:);
level.generationMeta = generationMeta;
level.samples2D = zeros(size(X, 1), 2);
level.projectionLabel = "";
end

function subsetLevels = buildDiagnosticArtifact(levelRecords, thresholds, Pf, beta, CoV, opts, varNames)
subsetLevels = struct();
subsetLevels.version = '1.0';
subsetLevels.createdAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
subsetLevels.failureDefinition = 'g(x) <= 0';
subsetLevels.Pf = Pf;
subsetLevels.beta = beta;
subsetLevels.CoV = CoV;
subsetLevels.N = opts.N;
subsetLevels.p0 = opts.p0;
subsetLevels.maxLevels = opts.maxLevels;
subsetLevels.varNames = varNames;
subsetLevels.levelThresholds = thresholds(:);
subsetLevels.thresholds = thresholds(:);
subsetLevels.levelCount = numel(levelRecords);
subsetLevels.levels = levelRecords;
subsetLevels.outDir = opts.outDir;
subsetLevels.space = 'raw_input_space';
subsetLevels.plotBasis = struct();
subsetLevels.metadata = struct();
subsetLevels.metadata.recommendedVariablePair = recommendVariablePair(levelRecords, numel(varNames));
subsetLevels.metadata.variableInfluence = estimateVariableInfluence(levelRecords, numel(varNames));
subsetLevels.metadata.description = strjoin({ ...
    'Each SS-k cloud represents the samples retained/generated at level k.', ...
    'The intermediate threshold b_k separates the subset event g(x) <= b_k.', ...
    'The final contour g(x)=0 is the failure boundary used in the reliability estimate.'}, ' ');
subsetLevels.metadata.failureBoundaryReconstruction = ...
    'The diagnostic figure reconstructs the displayed g(x)=0 contour in 2D with scatteredInterpolant using the saved samples and g-values.';
end

function idxPair = recommendVariablePair(levelRecords, nVars)
if nVars == 1
    idxPair = [1 1];
    return;
end
scores = estimateVariableInfluence(levelRecords, nVars);
[~, order] = sort(scores, 'descend');
idxPair = sort(order(1:min(2, nVars)));
if numel(idxPair) < 2
    idxPair = [idxPair(1) min(nVars, idxPair(1) + 1)];
end
end

function scores = estimateVariableInfluence(levelRecords, nVars)
allX = [];
allG = [];
for i = 1:numel(levelRecords)
    allX = [allX; levelRecords(i).samplesX]; %#ok<AGROW>
    allG = [allG; levelRecords(i).g(:)]; %#ok<AGROW>
end

scores = zeros(1, nVars);
for j = 1:nVars
    xj = allX(:, j);
    mask = isfinite(xj) & isfinite(allG);
    if nnz(mask) >= 3 && std(xj(mask)) > 0 && std(allG(mask)) > 0
        c = corrcoef(xj(mask), allG(mask));
        scores(j) = abs(c(1, 2));
    end
end
end

function subsetLevels = attachStored2DBasis(subsetLevels, opts)
basis = resolvePlotBasis(subsetLevels, opts);

for i = 1:subsetLevels.levelCount
    Xi = subsetLevels.levels(i).samplesX;
    subsetLevels.levels(i).samples2D = projectSamples(Xi, basis);
    subsetLevels.levels(i).projectionLabel = basis.label;
end

subsetLevels.plotBasis = basis;
end

function basis = resolvePlotBasis(subsetLevels, opts)
method = 'auto';
if isfield(opts, 'plotBasis') && isfield(opts.plotBasis, 'method') && ~isempty(opts.plotBasis.method)
    method = lower(string(opts.plotBasis.method));
end

allX = [];
for i = 1:subsetLevels.levelCount
    allX = [allX; subsetLevels.levels(i).samplesX]; %#ok<AGROW>
end
nVars = size(allX, 2);

if nVars == 1
    basis = struct('method', 'pair', 'variablePair', [1 1], ...
        'varNames', subsetLevels.varNames([1 1]), ...
        'label', sprintf('%s vs %s', subsetLevels.varNames(1), subsetLevels.varNames(1)));
    return;
end

if method == "pair" || (method == "auto" && hasVariablePair(opts))
    variablePair = parseVariablePair(opts.plotBasis.variablePair, subsetLevels.varNames);
    basis = struct('method', 'pair', 'variablePair', variablePair, ...
        'varNames', subsetLevels.varNames(variablePair), ...
        'label', sprintf('%s vs %s', subsetLevels.varNames(variablePair(1)), subsetLevels.varNames(variablePair(2))));
    return;
end

if method == "auto"
    variablePair = subsetLevels.metadata.recommendedVariablePair;
    basis = struct('method', 'pair', 'variablePair', variablePair, ...
        'varNames', subsetLevels.varNames(variablePair), ...
        'label', sprintf('%s vs %s', subsetLevels.varNames(variablePair(1)), subsetLevels.varNames(variablePair(2))));
    return;
end

[coeff, mu, sigma] = computePcaBasis(allX);
basis = struct();
basis.method = 'pca';
basis.coeff = coeff;
basis.mu = mu;
basis.sigma = sigma;
basis.label = 'PCA projection (PC1 vs PC2)';
basis.varNames = ["PC1"; "PC2"];
end

function tf = hasVariablePair(opts)
tf = isfield(opts, 'plotBasis') && isfield(opts.plotBasis, 'variablePair') ...
    && ~isempty(opts.plotBasis.variablePair);
end

function pair = parseVariablePair(variablePair, varNames)
if isnumeric(variablePair)
    pair = reshape(variablePair, 1, []);
else
    pair = zeros(1, numel(variablePair));
    for i = 1:numel(variablePair)
        pair(i) = find(string(varNames) == string(variablePair{i}), 1, 'first');
    end
end
if numel(pair) ~= 2 || any(isnan(pair)) || any(pair < 1) || any(pair > numel(varNames))
    error('A9_subset_simulation_pf:InvalidVariablePair', ...
        'plotBasis.variablePair must identify exactly two valid variables.');
end
pair = double(pair);
end

function [coeff, mu, sigma] = computePcaBasis(X)
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

function Z = projectSamples(X, basis)
switch lower(char(basis.method))
    case 'pair'
        if basis.variablePair(1) == basis.variablePair(2)
            Z = [X(:, basis.variablePair(1)), zeros(size(X, 1), 1)];
        else
            Z = X(:, basis.variablePair);
        end
    case 'pca'
        Xs = (X - basis.mu) ./ basis.sigma;
        Z = Xs * basis.coeff;
        if size(Z, 2) < 2
            Z(:, 2) = 0;
        end
    otherwise
        error('A9_subset_simulation_pf:UnknownBasis', 'Unsupported plot basis method.');
end
end

function covPf = subsetSimulationCoV(p0, N, nLevels, pfLast)
pfLast = max(min(pfLast, 1 - eps), eps);
termIntermediate = max(nLevels - 1, 0) * (1 - p0) / max(N * p0, eps);
termLast = (1 - pfLast) / max(N * pfLast, eps);
covPf = sqrt(max(termIntermediate + termLast, 0));
end

function beta = pf2beta(pf)
epsv = 1e-15;
pf = min(max(pf, epsv), 1 - epsv);
beta = -sqrt(2) * erfcinv(2 * pf);
end

function outSS = A9_run_subset_simulation_surrogate(opts)
% A9_RUN_SUBSET_SIMULATION_SURROGATE
% Driver para Subset Simulation usando o surrogate do stage 3:
%   g(X) = recalque_lim - y_hat(X)
%
% O script procura o fluxo em ./out_incremental quando executado a partir
% do projeto e tambem funciona em C:\Kazuo-Script.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end

    clc;

    work_dir = local_resolve_work_dir();
    out_dir = fullfile(work_dir, 'out_incremental');
    stage3_mat = fullfile(out_dir, 'stage3_best.mat');
    result_mat = fullfile(out_dir, 'A9_subset_simulation_result.mat');

    assert(isfile(stage3_mat), 'Nao encontrado: %s', stage3_mat);
    if exist('uq_evalModel', 'file') ~= 2
        error('A9_run_subset_simulation_surrogate:UQLabMissing', ...
            'uq_evalModel nao encontrado. Inicialize o UQLab antes de rodar o A9.');
    end
    if exist('uqlab', 'file') == 2
        try
            uqlab('-nosplash');
        catch
        end
    end

    S = load(stage3_mat, 'bestModel', 'muY', 'sdY', 'vars_train', 'Xtr');
    requiredFields = {'bestModel', 'muY', 'sdY', 'Xtr'};
    for i = 1:numel(requiredFields)
        if ~isfield(S, requiredFields{i})
            error('A9_run_subset_simulation_surrogate:MissingField', ...
                'Campo ausente em stage3_best.mat: %s', requiredFields{i});
        end
    end

    mdl = S.bestModel;
    muY = S.muY;
    sdY = S.sdY;
    Xtr = S.Xtr;

    if isfield(S, 'vars_train') && ~isempty(S.vars_train)
        vars_train = S.vars_train;
    else
        vars_train = arrayfun(@(k) sprintf('X%d', k), 1:size(Xtr, 2), 'UniformOutput', false);
    end

    if ~isfield(opts, 'recalque_lim') || isempty(opts.recalque_lim)
        opts.recalque_lim = 0.1;
    end
    if ~isfield(opts, 'N') || isempty(opts.N), opts.N = 1000; end
    if ~isfield(opts, 'p0') || isempty(opts.p0), opts.p0 = 0.10; end
    if ~isfield(opts, 'maxLevels') || isempty(opts.maxLevels), opts.maxLevels = 10; end
    if ~isfield(opts, 'proposalScale') || isempty(opts.proposalScale), opts.proposalScale = 0.8; end
    if ~isfield(opts, 'seed') || isempty(opts.seed), opts.seed = 123; end

    mu = mean(Xtr, 1, 'omitnan');
    sg = std(Xtr, 0, 1, 'omitnan');
    sg(~isfinite(sg) | sg <= 0) = 1e-6;

    gfunX = @(X) local_gfun_surrogate(X, mdl, muY, sdY, opts.recalque_lim);

    M = size(Xtr, 2);
    myInput = struct();
    myInput.Marginals = repmat(struct('Type', 'Gaussian', 'Parameters', [0 1]), 1, M);
    for k = 1:M
        myInput.Marginals(k).Type = 'Gaussian';
        myInput.Marginals(k).Parameters = [mu(k), sg(k)];
    end

    outSS = A9_subset_simulation_pf(gfunX, myInput, opts);
    outSS.recalque_lim = opts.recalque_lim;
    outSS.work_dir = work_dir;
    outSS.out_dir = out_dir;
    outSS.vars_train = vars_train;
    outSS.muX = mu;
    outSS.sigmaX = sg;

    save(result_mat, 'outSS', 'opts', 'mu', 'sg', 'vars_train');

    fprintf('\n=== Subset Simulation via surrogate ===\n');
    fprintf('Diretorio de trabalho: %s\n', work_dir);
    fprintf('Pf   = %.8g\n', outSS.Pf);
    fprintf('beta = %.6f\n', outSS.beta);
    fprintf('CoV  = %.6f\n', outSS.CoV);
    fprintf('nLv  = %d\n', outSS.nLevels);
    fprintf('failReached = %d\n', outSS.failReached);
    fprintf('Resultado salvo em: %s\n', result_mat);
end

function g = local_gfun_surrogate(X, mdl, muY, sdY, recalque_lim)
    yN_hat = uq_evalModel(mdl, X);
    y_hat = yN_hat * sdY + muY;
    g = recalque_lim - y_hat;
    g = g(:);
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
        if isfolder(fullfile(base, 'out_incremental')) || isfile(fullfile(base, 'out_incremental', 'stage3_best.mat'))
            work_dir = base;
            return;
        end
    end

    error('A9_run_subset_simulation_surrogate:WorkDirNotFound', ...
        'Nao foi possivel localizar o diretorio do workflow (out_incremental).');
end

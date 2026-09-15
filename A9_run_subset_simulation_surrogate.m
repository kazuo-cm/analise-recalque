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

    if ~isfile(stage3_mat)
        error('A9_run_subset_simulation_surrogate:MissingStage3', ...
            'Nao encontrado: %s', stage3_mat);
    end
    if exist('uqlab', 'file') == 2
        try
            uqlab('-nosplash');
        catch
        end
    end
    if exist('uq_evalModel', 'file') ~= 2
        error('A9_run_subset_simulation_surrogate:UQLabMissing', ...
            'uq_evalModel nao encontrado. Inicialize o UQLab antes de rodar o A9.');
    end

    S = load(stage3_mat);
    requiredFields = {'bestModel', 'muY', 'sdY'};
    for i = 1:numel(requiredFields)
        if ~isfield(S, requiredFields{i})
            error('A9_run_subset_simulation_surrogate:MissingField', ...
                'Campo ausente em stage3_best.mat: %s', requiredFields{i});
        end
    end

    mdl = S.bestModel;
    muY = S.muY;
    sdY = S.sdY;
    myInput = local_resolve_input_model(S, opts);
    M = numel(myInput.Marginals);

    if isfield(S, 'vars_train') && ~isempty(S.vars_train)
        vars_train = S.vars_train;
    else
        vars_train = arrayfun(@(k) sprintf('X%d', k), 1:M, 'UniformOutput', false);
    end

    if isfield(S, 'Xtr') && ~isempty(S.Xtr) && size(S.Xtr, 2) ~= M
        error('A9_run_subset_simulation_surrogate:DimensionMismatch', ...
            'Xtr possui %d colunas, mas myInput tem %d marginais.', size(S.Xtr, 2), M);
    end

    if ~isfield(opts, 'recalque_lim') || isempty(opts.recalque_lim)
        opts.recalque_lim = 0.1;
    end
    if ~isfield(opts, 'N') || isempty(opts.N), opts.N = 1000; end
    if ~isfield(opts, 'p0') || isempty(opts.p0), opts.p0 = 0.10; end
    if ~isfield(opts, 'maxLevels') || isempty(opts.maxLevels), opts.maxLevels = 10; end
    if ~isfield(opts, 'proposalScale') || isempty(opts.proposalScale), opts.proposalScale = 0.8; end
    if ~isfield(opts, 'seed') || isempty(opts.seed), opts.seed = 123; end
    if ~isfield(opts, 'recordDiagnostics') || isempty(opts.recordDiagnostics), opts.recordDiagnostics = true; end

    gfunX = @(X) local_gfun_surrogate(X, mdl, muY, sdY, opts.recalque_lim);

    outSS = A9_subset_simulation_pf(gfunX, myInput, opts);
    outSS.recalque_lim = opts.recalque_lim;
    outSS.work_dir = work_dir;
    outSS.out_dir = out_dir;
    outSS.vars_train = vars_train;
    outSS.myInput = myInput;

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    save(result_mat, 'outSS', 'opts', 'myInput', 'vars_train');

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
    y_hat = bsxfun(@plus, bsxfun(@times, yN_hat, sdY), muY);
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

function myInput = local_resolve_input_model(S, opts)
    if isfield(opts, 'myInput') && local_has_marginals(opts.myInput)
        myInput = opts.myInput;
        return;
    end

    stageCandidates = {'myInput', 'uqInput', 'inputModel'};
    for i = 1:numel(stageCandidates)
        name = stageCandidates{i};
        if isfield(S, name) && local_has_marginals(S.(name))
            myInput = S.(name);
            return;
        end
    end

    if isfield(opts, 'muX') && isfield(opts, 'sigmaX')
        mu = opts.muX(:)';
        sg = opts.sigmaX(:)';
        if numel(mu) ~= numel(sg)
            error('A9_run_subset_simulation_surrogate:InvalidMoments', ...
                'opts.muX e opts.sigmaX devem ter o mesmo numero de elementos.');
        end
        if any(~isfinite(sg) | sg <= 0)
            error('A9_run_subset_simulation_surrogate:InvalidMoments', ...
                'opts.sigmaX contem valores invalidos (todos devem ser finitos e > 0).');
        end
        myInput = struct();
        myInput.Marginals = repmat(struct('Type', 'Gaussian', 'Parameters', [0 1]), 1, numel(mu));
        for k = 1:numel(mu)
            myInput.Marginals(k).Type = 'Gaussian';
            myInput.Marginals(k).Parameters = [mu(k), sg(k)];
        end
        return;
    end

    error('A9_run_subset_simulation_surrogate:MissingInputModel', ...
        ['Nao foi possivel reconstruir as marginais do problema. ' ...
         'Forneca opts.myInput (preferencialmente) ou opts.muX/opts.sigmaX.']);
end

function tf = local_has_marginals(candidate)
    tf = isstruct(candidate) && isfield(candidate, 'Marginals') && ~isempty(candidate.Marginals);
end

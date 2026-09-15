function outSS = A9_subset_simulation_pf(gfunX, myInput, opts)
% A9_SUBSET_SIMULATION_PF
% Estima Pf = P[g(X)<=0] via Subset Simulation (Au & Beck, 2001).

    if nargin < 2
        error('A9_subset_simulation_pf:NotEnoughInputs', ...
            'Forneca gfunX e myInput.');
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'N') || isempty(opts.N), opts.N = 1000; end
    if ~isfield(opts, 'p0') || isempty(opts.p0), opts.p0 = 0.10; end
    if ~isfield(opts, 'maxLevels') || isempty(opts.maxLevels), opts.maxLevels = 10; end
    if ~isfield(opts, 'proposalScale') || isempty(opts.proposalScale), opts.proposalScale = 0.8; end
    if isfield(opts, 'seed') && ~isempty(opts.seed)
        rng(opts.seed, 'twister');
    end

    validateattributes(gfunX, {'function_handle'}, {'scalar'}, mfilename, 'gfunX', 1);
    if ~isstruct(myInput) || ~isfield(myInput, 'Marginals') || isempty(myInput.Marginals)
        error('A9_subset_simulation_pf:InvalidInput', ...
            'myInput.Marginals deve conter as marginais Gaussianas.');
    end

    N = opts.N;
    p0 = opts.p0;
    Lmax = opts.maxLevels;
    sProp = opts.proposalScale;

    validateattributes(N, {'numeric'}, {'scalar', 'integer', '>=', 10}, mfilename, 'opts.N');
    validateattributes(p0, {'numeric'}, {'scalar', '>', 0, '<', 1}, mfilename, 'opts.p0');
    validateattributes(Lmax, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'opts.maxLevels');
    validateattributes(sProp, {'numeric'}, {'scalar', 'positive'}, mfilename, 'opts.proposalScale');

    nKeepFloat = p0 * N;
    nKeep = round(nKeepFloat);
    if abs(nKeepFloat - nKeep) > 1e-10
        error('A9_subset_simulation_pf:InvalidSubsetFraction', ...
            'Escolha opts.p0 e opts.N tais que p0*N seja inteiro. Recebido: %.12g.', nKeepFloat);
    end
    nKeep = max(1, nKeep);
    p0 = nKeep / N;

    M = numel(myInput.Marginals);
    u2x = @(U) local_u2x_gauss(U, myInput);

    U = randn(N, M);
    X = u2x(U);
    g = gfunX(X);
    g = g(:);

    if numel(g) ~= N
        error('A9_subset_simulation_pf:InvalidResponse', ...
            'gfunX deve retornar um valor por linha de X.');
    end

    bLevels = nan(Lmax, 1);
    level = 0;
    failReached = false;

    while level < Lmax
        level = level + 1;

        gs = sort(g, 'ascend');
        b = gs(nKeep);
        bLevels(level) = b;

        if b <= 0
            failReached = true;
            break;
        end

        idxSeed = find(g <= b);
        if numel(idxSeed) < nKeep
            [~, ord] = sort(g, 'ascend');
            idxSeed = ord(1:nKeep);
        else
            idxSeed = idxSeed(1:nKeep);
        end

        Useed = U(idxSeed, :);
        gseed = g(idxSeed);

        Unew = zeros(N, M);
        gnew = zeros(N, 1);

        chainLen = ceil(N / nKeep);
        c = 0;
        for i = 1:nKeep
            uc = Useed(i, :);
            gc = gseed(i);

            for t = 1:chainLen
                up = uc + sProp * randn(1, M);
                gp = gfunX(u2x(up));
                gp = gp(:);
                if numel(gp) ~= 1
                    error('A9_subset_simulation_pf:InvalidResponse', ...
                        'gfunX deve retornar um unico valor para cada amostra proposta.');
                end
                gp = gp(1);

                if gp <= b
                    uc = up;
                    gc = gp;
                end

                c = c + 1;
                if c <= N
                    Unew(c, :) = uc;
                    gnew(c, 1) = gc;
                else
                    break;
                end
            end

            if c >= N
                break;
            end
        end

        if c < N
            error('A9_subset_simulation_pf:ChainFillFailed', ...
                'Numero insuficiente de amostras condicionais geradas no nivel %d.', level);
        end

        U = Unew(1:N, :);
        g = gnew(1:N);
    end

    nLevels = level;
    pLast = mean(g <= 0);
    nIntermediate = max(nLevels - 1, 0);
    if ~failReached
        nIntermediate = nLevels;
    end

    Pf = (p0 ^ nIntermediate) * pLast;
    Pf = max(min(Pf, 1 - 1e-15), 1e-15);
    beta = local_beta_from_pf(Pf);

    if pLast > 0
        relVarLast = max(1 - pLast, 0) / (N * pLast);
    else
        relVarLast = NaN;
    end

    if nIntermediate == 0
        CoV = sqrt(relVarLast);
    else
        relVarSubset = nIntermediate * (1 - p0) / (N * p0);
        if isnan(relVarLast)
            CoV = NaN;
        else
            CoV = sqrt(relVarSubset + relVarLast);
        end
    end

    outSS = struct();
    outSS.Pf = Pf;
    outSS.beta = beta;
    outSS.CoV = CoV;
    outSS.nLevels = nLevels;
    outSS.bLevels = bLevels(1:nLevels);
    outSS.p0 = p0;
    outSS.N = N;
    outSS.failReached = failReached;
    outSS.pLast = pLast;
end

function X = local_u2x_gauss(U, myInput)
    [N, M] = size(U);
    if M ~= numel(myInput.Marginals)
        error('A9_subset_simulation_pf:DimensionMismatch', ...
            'U possui %d colunas, mas myInput contem %d marginais.', M, numel(myInput.Marginals));
    end
    X = zeros(N, M);
    for k = 1:M
        if isfield(myInput.Marginals(k), 'Type') && ~isempty(myInput.Marginals(k).Type)
            thisType = lower(string(myInput.Marginals(k).Type));
            if ~(thisType == "gaussian" || thisType == "normal")
                error('A9_subset_simulation_pf:UnsupportedMarginal', ...
                    'Marginal %d do tipo %s nao eh suportada por esta implementacao gaussiana.', ...
                    k, char(thisType));
            end
        end
        params = myInput.Marginals(k).Parameters;
        if numel(params) < 2
            error('A9_subset_simulation_pf:InvalidMarginal', ...
                'Marginal %d deve ter [mu sigma].', k);
        end
        mu = params(1);
        sg = params(2);
        X(:, k) = mu + sg * U(:, k);
    end
end

function beta = local_beta_from_pf(Pf)
    Pf = max(min(Pf, 1 - 1e-15), 1e-15);
    beta = -sqrt(2) * erfcinv(2 * Pf);
end

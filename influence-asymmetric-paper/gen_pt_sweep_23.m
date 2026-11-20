%% Instructions

% Need to do before running:
% Set "desiredWorkers" in parallel section based on hardware
% Set interval and number of gen pts to sweep over

%% Set parameters

% Number of gen points
n = 750;

% Set fixed number of clusters expected
kclus = 3; % <-- set expected number of clusters here

% Model parameters, Signaling and Responsive Regions
eps = .01;              % Introduce gap between bud_pt and s1
bud_pt = 0.6;           % position at which budding takes place
r1 = bud_pt - (1 - .8);   % beginning of responsive region
r2 = bud_pt;                 % end of responsive region
s1 = bud_pt + eps;           % beginning of S phase
s2 = s1 + .27;             % end of S Phase
beta = -0.6;                 % recruitment(+)/retardation(-) rate
gen_pt = linspace(.28,.4,n); % values of g to loop over

% Original model parameters
Cap = 10200;       % Carrying capacity
gamma = .51;        % coefficient of logistic nonlinearity

% Simulation parameters
ncells0 = Cap;      % initial number of cells
bins = 98;          % number of bins used for histogram
cycles = 300;       % number of cell cycles
h = .004;           % step size
steps = ceil(cycles/h);
sigma = .0000000005;      % noise level from main script
sqrh = sqrt(h);
hbeta = h*beta;

% Results arrays
MaxDensity = zeros(size(gen_pt));
MinDensity = zeros(size(gen_pt));

MaxBud = zeros(size(gen_pt));
MinBud = zeros(size(gen_pt));

ClusterSD = zeros(size(gen_pt));

% Average generation fractions over the last 2 cycles
AvgD  = zeros(size(gen_pt));
AvgM1 = zeros(size(gen_pt));
AvgM2 = zeros(size(gen_pt));
AvgM3 = zeros(size(gen_pt));
AvgM4 = zeros(size(gen_pt));
AvgM5 = zeros(size(gen_pt));
AvgM6 = zeros(size(gen_pt));

%% Parallel stuff
maxNumCompThreads(1);    % disable internal BLAS threading
desiredWorkers = 24;      % CHECK and set

% Update cluster max workers if needed
c = parcluster('local');
if c.NumWorkers < desiredWorkers
    c.NumWorkers = desiredWorkers;
    saveProfile(c);
end

% Restart pool if needed
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= desiredWorkers
    delete(gcp('nocreate'));
    parpool('Processes', desiredWorkers);
end

fprintf('Parallel pool started with %d workers.\n', gcp().NumWorkers);

% Start timing
overall_tic = tic;

%% Compute
parfor k = 1:length(gen_pt)

    % Use the same harvest-rate relation as in your previous g-loop script
    % so the loop still scans the same effective parameterization.
    mu_k = log((1+sqrt(5))/2)/(1-gen_pt(k));
    harvest_step = (1 - exp(-mu_k*h));

    % Per-worker RNG
    rng(k,'twister');

    % Initial condition: same style as main script
    u = unifrnd(0,1,ncells0,1);

    % Generation labels
    ug = zeros(size(u));
    ug(2:2:end) = 1;
    ug(4:4:end) = 2;
    ug(8:8:end) = 3;
    ug(16:16:end) = 4;

    % Time/history variables
    t = 0;
    BudIdx = [];
    T = [];
    TotCells = [];

    % Histories of generation fractions over the final 2 cycles
    D_hist  = [];
    M1_hist = [];
    M2_hist = [];
    M3_hist = [];
    M4_hist = [];
    M5_hist = [];
    M6_hist = [];

    for i = 1:steps
        ncells = length(u);

        % Logistic growth term from main script
        u = u + h*(1 + gamma*(Cap - ncells)/Cap);
        t = t + h;

        % Cells in S phase
        S = find(s1 <= u & u <= s2);
        I = numel(S)/ncells;

        % Cells in responsive region
        J = find(r1 <= u & u <= r2);
        u(J) = u(J) + hbeta*I*ones(size(J));

        % Add noise
        u = u + sqrh*sigma*(rand(1,ncells)-.5)';

        % Division
        Div = find(u > 1);
        if ~isempty(Div)
            u_d = u(Div) - 1;                  % positions of new daughters
            u(Div) = u(Div) + gen_pt(k) - 1;  % updated mothers
            ug(Div) = ug(Div) + 1;            % mothers increment generation

            % Append daughters
            u = [u; u_d];
            ug = [ug; zeros(size(Div))];
        end

        % Harvest
        keep = rand(size(u)) > harvest_step;
        u = u(keep);
        ug = ug(keep);

        ncells = length(u);
        totcells = ncells;

        % Record during the final 2 cycles
        if t >= cycles - 2
            TotCells = [TotCells, totcells];

            Bud = find(u > bud_pt);
            bi = length(Bud)/totcells;
            BudIdx = [BudIdx; bi];
            T = [T; t];

            d  = sum(ug == 0)/totcells;
            m1 = sum(ug == 1)/totcells;
            m2 = sum(ug == 2)/totcells;
            m3 = sum(ug == 3)/totcells;
            m4 = sum(ug == 4)/totcells;
            m5 = sum(ug == 5)/totcells;
            m6 = sum(ug >= 6)/totcells;

            D_hist  = [D_hist; d];
            M1_hist = [M1_hist; m1];
            M2_hist = [M2_hist; m2];
            M3_hist = [M3_hist; m3];
            M4_hist = [M4_hist; m4];
            M5_hist = [M5_hist; m5];
            M6_hist = [M6_hist; m6];
        end
    end

    % Save density and bud summaries
    MaxDensity(k) = max(TotCells);
    MinDensity(k) = min(TotCells);

    MaxBud(k) = max(BudIdx);
    MinBud(k) = min(BudIdx);

    % Save average generation fractions over final 2 cycles
    if ~isempty(D_hist)
        AvgD(k)  = mean(D_hist);
        AvgM1(k) = mean(M1_hist);
        AvgM2(k) = mean(M2_hist);
        AvgM3(k) = mean(M3_hist);
        AvgM4(k) = mean(M4_hist);
        AvgM5(k) = mean(M5_hist);
        AvgM6(k) = mean(M6_hist);
    end

    %% Compute final cluster spread statistic
    % Use final cell positions, mapped to unit circle
    u_col = mod(u(:),1);
    thetas = u_col * 2*pi;
    X = [cos(thetas), sin(thetas)];

    % Cluster using fixed number of clusters
    if size(X,1) <= kclus
        idx_h = (1:size(X,1))';
    else
        opts = statset('MaxIter',500,'Display','off');
        idx_h = kmeans(X, kclus, 'Replicates', 5, 'Options', opts);
    end

    % Circular-aware standard deviation within each cluster
    clusStd = zeros(kclus,1);
    clusCount = zeros(kclus,1);

    for cidx = 1:kclus
        mask = (idx_h == cidx);
        vals = u_col(mask);
        clusCount(cidx) = numel(vals);

        if clusCount(cidx) <= 1
            clusStd(cidx) = 0;
        else
            thetas_c = vals * 2*pi;
            C = mean(cos(thetas_c));
            S = mean(sin(thetas_c));
            R = sqrt(C.^2 + S.^2);
            circStd = sqrt(-2*log(max(R,eps))); % radians
            clusStd(cidx) = circStd / (2*pi);   % convert to cycle fraction
        end
    end

    totalCellsFinal = sum(clusCount);
    if totalCellsFinal > 0
        weightedAvgStd = sum(clusStd .* (clusCount/totalCellsFinal));
    else
        weightedAvgStd = 0;
    end

    ClusterSD(k) = weightedAvgStd;
end

%% End timing
elapsed_time = toc(overall_tic);
fprintf('Done! Total runtime: %.2f seconds (%.2f minutes)\n', elapsed_time, elapsed_time/60);

%% Save workspace as .mat file
timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');
filename = ['gloop_cellwise_wsp_' timestamp '.mat'];
save(filename)
%% MCMC parameter optimization for gas-pocket forward model (two samples)
% -------------------------------------------------------------------------
% Purpose
%   Estimate physical/source parameters by fitting modeled pressure spectra
%   to observed spectra using Metropolis–Hastings MCMC.
%
% Data / forward model
%   - Observations: two amplitude spectra (sample #1 and #2), stored in
%       spect_store_100(:,1) and spect_store_100(:,2)
%   - Forward model: forwardmodelgas15_10ps_ex2(m, t0, qn, sigma)
%       returns:
%         Posc_in  : time-domain pressure (unused here except for length)
%         Pspect   : complex pressure spectrum
%
% Two-sample setup (shared / coupled prior)
%   - Sample #1: "9 capillary" (parameter vector m)
%   - Sample #2: "dense sponge" (parameter vector m2)
%   - Each sample has its own (kappa, phi, Q, t0, sigma, qn_n)
%   - Coupling: sigma_diff = sigma - sigma_d is penalized in the prior
%
% Key parameters
%   kk        : run label (controls output filenames)
%   tau (s)   : signal duration
%   dt (s)    : sampling interval for forward-model time series
%   fc (Hz)   : sensor sampling frequency (informational, not required below)
%
% Units / conversions
%   Q is in kg/s of air; conversion from ml/s is:
%     Q[kg/s] = (ml/s) * 1.2e-6
%     (1.2e-6 kg/s per ml/s is used here)
%
% Gas excitation parameterization
%   - NN impulses (fixed)
%   - t0 : impulse times (seconds), length NN
%   - qn_n : impulse weights (dimensionless ratios), length NN
%   - qn = qn_n * Q * tau (kg) gives the amount of gas for each impulse
%
% Priors (examples; adjust to your study)
%   - phi (sample1): Normal(mean=0.006, sd=0.01)
%   - phi (sample2): Normal(mean=0.9,   sd=0.1)
%   - sigma coupling: Normal(0, sd=0.0005) on (sigma - sigma_d)
%   - qn_n: Gaussian penalty around its mean to keep weights near-uniform
%
% Proposal strategy
%   A simple component-wise schedule using mod(k,8) to perturb one block
%   of variables at a time:
%     0: (kappa, phi, Q) for sample1
%     1: sigma for sample1
%     2: t0 for sample1
%     3: (kappa, phi, Q) for sample2
%     4: t0 for sample2
%     5: sigma for sample2
%     6: qn_n for sample1
%     7: qn_n for sample2
%
% Dependencies (must be on MATLAB path)
%   - forwardmodelgas15_10ps_ex2
%   - edgeAdjustedAverageFilter2
%   - generategas_3, generategas_6
%   - generate_sig2, generateexc_3, generateqn
%


clear; clc;

%% --------------------------- Paths & data -------------------------------
% (Keep paths minimal in publishable code; prefer relative paths or add
%  instructions in README.)
addpath('/Users/kyuninkim/Downloads');
addpath('/Users/kyuninkim/Downloads/Pdata_inv_joint');
addpath('/Users/kyuninkim/Downloads/GA_parameter_2');
addpath('/Users/kyuninkim/Downloads/Numerical_Simu');

addpath('/import/c1/SOLIDEARTH/kkim39/exp_inv_store_50');
addpath('/import/c1/SOLIDEARTH/kkim39/exp_inv_store_100');
addpath('/import/c1/SOLIDEARTH/kkim39/Lib_exp');

load('Spect_store_50ml_new.mat','f_h');                % frequency vector
load('spect_store_100_store2.mat','spect_store_100');  % observed spectra

% Observed amplitude spectra (row vectors)
obs_amp_1 = spect_store_100(:,1)';   % sample #1
obs_amp_2 = spect_store_100(:,2)';   % sample #2

%% ---------------------- MCMC configuration ------------------------------
kk  = 22;        % run label for output files
tau = 10;        % signal duration (s)
niter = 2_000_000;

% Forward-model time discretization
dt = 0.0002;              % (s)
t  = 0:dt:tau; %#ok<NASGU> % (s) for reference

% Sensor sampling frequency (informational)
fc = 5000; %#ok<NASGU>

% --- Global step sizes used by your proposal functions ---
% Prefer a struct in publishable code, but keep global if your proposal
% functions depend on it.
global Stepsize1 Stepsize2 Stepsize3 Stepsize4 Stepsize5 Stepsize6
Stepsize1 = 5e-5;
Stepsize2 = 2e-3;
Stepsize3 = 1.5e-7;
Stepsize4 = 2e-4;
Stepsize5 = 4e-7;
Stepsize6 = 9e-4;

%% ------------------------ Fixed parameters ------------------------------
T     = 293.15;   % temperature (K)
R     = 0.02;     % conduit radius (m)
D     = 0.4;      % gas pocket thickness/length scale (m) (fixed here)
L_1   = 0.10;     % sample #1 porous length (m)  (was L_val2)
L_2   = 0.20;     % sample #2 porous length (m)  (was L_val)

% NOTE: remove unused parameters (rho_s, Qf) as requested

%% ------------------ Initial state (two samples) -------------------------
NN = 20;                 % number of impulses (fixed)
sigma0   = 0.013;        % initial impulse timescale for sample #1 (s)
sigma0_d = 0.013;        % initial impulse timescale for sample #2 (s)

% Initial permeability parameters in log10 space
kappa0   = -8.1;         % log10(permeability) for sample #1
kappa0_d = -8.1;         % log10(permeability) for sample #2

% Initial porosities (dimensionless)
phi0   = 0.06;           % sample #1
phi0_d = 0.45;           % sample #2

% Initial gas flow rates (kg/s). 1.2e-6 converts ml/s -> kg/s (air).
Q_ml_s = 80;             % initial flow in ml/s (for readability)
Q0     = Q_ml_s * 1.2e-6;

% Initial impulse times (uniform in [0, tau])
t0   = rand(1,NN) * tau;     % sample #1
t0_d = rand(1,NN) * tau;     % sample #2

% Initial impulse weights qn_n (dimensionless ratios that sum ~ 1)
% Here: 20% standard deviation around uniform, mean-corrected
qn_std_frac = 0.20;
qn_var   = randn(1,NN) * (1/NN) * qn_std_frac;  qn_var   = qn_var - mean(qn_var);
qn_var_d = randn(1,NN) * (1/NN) * qn_std_frac;  qn_var_d = qn_var_d - mean(qn_var_d);

qn_n   = (1/NN) * ones(1,NN) + qn_var;    % sample #1
qn_n_d = (1/NN) * ones(1,NN) + qn_var_d;  % sample #2

% Convert ratios -> gas mass per impulse (kg)
qn   = qn_n   * Q0 * tau;
qn_d = qn_n_d * Q0 * tau;

%% ---------------------- Prior definitions -------------------------------
% Prior on qn_n to keep weights near-uniform (Gaussian penalty around mean)
% NOTE: you described this as “qn should be gaussian distribution”.
qn_n_mean   = mean(qn_n);
qn_n_mean_d = mean(qn_n_d);

% Prior scale (keep your original numbers, but define clearly)
% In original code: denom = 2*(0.05*0.2)^2
qn_prior_sd = (0.05 * 0.2);

lp_qn   = @(x) -sum((x - mean(x)).^2) / (2 * qn_prior_sd^2);

% Priors on porosity (as in your original)
lp_phi_1 = @(phi) -(phi - 0.006)^2 / (2 * 0.01^2);
lp_phi_2 = @(phi) -(phi - 0.9  )^2 / (2 * 0.10^2);

% Coupling prior on sigma difference
lp_sigma_couple = @(s1,s2) -((s1 - s2)^2) / (2 * (0.0005^2));

%% ---------------------- Likelihood setup --------------------------------
% We model residuals as Gaussian with diagonal covariance.
% You used Std_signal = 0.5*obs_amp (elementwise), Cov = diag(Std_signal).
% For diagonal cov, inv(diag(v)) = diag(1./v). Use this for speed/stability.

idx = 10:801;  % frequency bins used in likelihood (as in original)

Std1 = 0.5 * obs_amp_1(:);              % column
Std2 = 0.5 * obs_amp_2(:);

Icov1 = diag(1 ./ Std1(idx));           % inverse diag covariance
Icov2 = diag(1 ./ Std2(idx));

% Smoothing settings for modeled amplitude spectrum
window_size = 31;

%% ------------------ Helper: run forward model & misfit ------------------
% Convert kappa (log10) -> permeability (m^2)
kappa_to_val = @(klog) 10^(klog);

% Pack parameter vector expected by forwardmodelgas15_10ps_ex2:
%   m = [T, Q, R, L, kappa, D, phi]
make_m = @(Q,R,L,kappa, D, phi) [T, Q, R, L, kappa, D, phi];

% Compute smoothed modeled amplitude spectrum
model_smooth_amp = @(Pspect) edgeAdjustedAverageFilter2(abs(Pspect), window_size);

% Log-likelihood for one sample given modeled smoothed spectrum
loglike_one = @(obs_amp, model_amp, Icov) ...
    (-1/2) * ( (obs_amp(idx) - model_amp(idx)) * Icov * (obs_amp(idx) - model_amp(idx))' );

%% ----------------------- Initialize current state ------------------------
cur.kappa   = kappa0;
cur.sigma   = sigma0;
cur.t0      = t0;
cur.Q       = Q0;
cur.phi     = phi0;
cur.qn_n    = qn_n;

cur2.kappa  = kappa0_d;
cur2.sigma  = sigma0_d;
cur2.t0     = t0_d;
cur2.Q      = Q0;
cur2.phi    = phi0_d;
cur2.qn_n   = qn_n_d;

% Current priors on qn_n
cur.lp_qn   = lp_qn(cur.qn_n);
cur2.lp_qn  = lp_qn(cur2.qn_n);

% Initial forward model evaluations
m1 = make_m(cur.Q,  R, L_1, kappa_to_val(cur.kappa),  D, cur.phi);
m2 = make_m(cur2.Q, R, L_2, kappa_to_val(cur2.kappa), D, cur2.phi);

[~, Pspect1, ~, ~] = forwardmodelgas15_10ps_ex2(m1, cur.t0,  cur.qn_n  * cur.Q  * tau, cur.sigma);
[~, Pspect2, ~, ~] = forwardmodelgas15_10ps_ex2(m2, cur2.t0, cur2.qn_n * cur2.Q * tau, cur2.sigma);

Pspect1(1) = 0;  Pspect2(1) = 0;

cur.model_amp  = model_smooth_amp(Pspect1);
cur2.model_amp = model_smooth_amp(Pspect2);

% Current log-likelihoods
ll1 = loglike_one(obs_amp_1, cur.model_amp,  Icov1);
ll2 = loglike_one(obs_amp_2, cur2.model_amp, Icov2);
ll_current = ll1 + ll2;

% Current log-prior (phi priors + sigma coupling + qn priors)
lp_current = lp_phi_1(cur.phi) + lp_phi_2(cur2.phi) + ...
             lp_sigma_couple(cur.sigma, cur2.sigma) + ...
             cur.lp_qn + cur2.lp_qn;

% MAP tracker (store maximum posterior = ll + lp)
post_best = -inf;

%% -------------------------- Storage arrays ------------------------------
kappastore   = zeros(1,niter);
sigmastore   = zeros(1,niter);
gasstore     = zeros(1,niter);
phistore     = zeros(1,niter);

kappastore2  = zeros(1,niter);
sigmastore2  = zeros(1,niter);
gasstore2    = zeros(1,niter);
phistore2    = zeros(1,niter);

llstore      = zeros(1,niter);
acceptance   = zeros(1,niter);

% Optional thinning storage (every 20000 steps -> 100 samples)
thin_step = 20000;
thinN = floor(niter/thin_step);
t0_store   = zeros(NN, thinN);
t0_store2  = zeros(NN, thinN);
qn_n_store  = zeros(NN, thinN);
qn_n_store2 = zeros(NN, thinN);

%% ============================= MCMC loop ================================
warning('off', 'MATLAB:colon:nonIntegerIndex');

for k = 1:niter

    % Decide which parameter block to perturb
    ka = rem(k,8);

    % Start from current state (candidate = current)
    cand  = cur;
    cand2 = cur2;

    % Propose a move (proposal functions are yours)
    if ka == 0
        % Sample #1: propose permeability/porosity/flow
        [cand.kappa, cand.phi, cand.Q] = generategas_3(cur.kappa, cur.phi, cur.Q);

    elseif ka == 1
        % Sample #1: propose sigma
        cand.sigma = generate_sig2(cur.sigma);

    elseif ka == 2
        % Sample #1: propose excitation times t0
        cand.t0 = generateexc_3(cur.t0);

    elseif ka == 3
        % Sample #2: propose permeability/porosity/flow
        [cand2.kappa, cand2.phi, cand2.Q] = generategas_6(cur2.kappa, cur2.phi, cur2.Q);

    elseif ka == 4
        % Sample #2: propose excitation times t0
        cand2.t0 = generateexc_3(cur2.t0);

    elseif ka == 5
        % Sample #2: propose sigma
        cand2.sigma = generate_sig2(cur2.sigma);

    elseif ka == 6
        % Sample #1: propose qn_n (weights) + its prior
        [cand.qn_n, cand.lp_qn] = generateqn(cur.qn_n);

    else
        % Sample #2: propose qn_n (weights) + its prior
        [cand2.qn_n, cand2.lp_qn] = generateqn(cur2.qn_n);
    end

    % Ensure priors for qn are defined if not updated in this move
    if ~isfield(cand,'lp_qn') || isempty(cand.lp_qn),   cand.lp_qn  = lp_qn(cand.qn_n);  end
    if ~isfield(cand2,'lp_qn') || isempty(cand2.lp_qn), cand2.lp_qn = lp_qn(cand2.qn_n); end

    % --- Forward model evaluation (only recompute what changed) ---
    % For simplicity and clarity (publishable), recompute the affected sample.
    % (You can optimize further by caching and reusing spectra.)

    ll_cand_1 = ll1;
    ll_cand_2 = ll2;
    cand_model_amp_1 = cur.model_amp;
    cand_model_amp_2 = cur2.model_amp;

    if ismember(ka, [0 1 2 6])
        % Update sample #1 forward model
        m1c = make_m(cand.Q, R, L_1, kappa_to_val(cand.kappa), D, cand.phi);
        qn_c = cand.qn_n * cand.Q * tau;

        [~, Ps1c, ~, ~] = forwardmodelgas15_10ps_ex2(m1c, cand.t0, qn_c, cand.sigma);
        Ps1c(1) = 0;

        cand_model_amp_1 = model_smooth_amp(Ps1c);
        ll_cand_1 = loglike_one(obs_amp_1, cand_model_amp_1, Icov1);
    end

    if ismember(ka, [3 4 5 7])
        % Update sample #2 forward model
        m2c = make_m(cand2.Q, R, L_2, kappa_to_val(cand2.kappa), D, cand2.phi);
        qn2_c = cand2.qn_n * cand2.Q * tau;

        [~, Ps2c, ~, ~] = forwardmodelgas15_10ps_ex2(m2c, cand2.t0, qn2_c, cand2.sigma);
        Ps2c(1) = 0;

        cand_model_amp_2 = model_smooth_amp(Ps2c);
        ll_cand_2 = loglike_one(obs_amp_2, cand_model_amp_2, Icov2);
    end

    ll_cand = ll_cand_1 + ll_cand_2;

    % Candidate log-prior
    lp_cand = lp_phi_1(cand.phi) + lp_phi_2(cand2.phi) + ...
              lp_sigma_couple(cand.sigma, cand2.sigma) + ...
              cand.lp_qn + cand2.lp_qn;

    % Metropolis–Hastings accept/reject
    logalpha = (ll_cand + lp_cand) - (ll_current + lp_current);
    logalpha = min(0, logalpha);

    if log(rand()) < logalpha
        % Accept
        cur  = cand;
        cur2 = cand2;

        ll1 = ll_cand_1;
        ll2 = ll_cand_2;

        ll_current = ll_cand;
        lp_current = lp_cand;

        cur.model_amp  = cand_model_amp_1;
        cur2.model_amp = cand_model_amp_2;

        acceptance(k) = 1;
    else
        % Reject
        acceptance(k) = 0;
    end

    % Store chain
    kappastore(k)  = cur.kappa;
    sigmastore(k)  = cur.sigma;
    gasstore(k)    = cur.Q;
    phistore(k)    = cur.phi;

    kappastore2(k) = cur2.kappa;
    sigmastore2(k) = cur2.sigma;
    gasstore2(k)   = cur2.Q;
    phistore2(k)   = cur2.phi;

    llstore(k)     = ll_current;

    % MAP (max posterior) tracking
    post_now = ll_current + lp_current;
    if post_now > post_best
        post_best = post_now;

        MAP.kappa  = cur.kappa;
        MAP.sigma  = cur.sigma;
        MAP.phi    = cur.phi;
        MAP.Q      = cur.Q;
        MAP.t0     = cur.t0;
        MAP.qn_n   = cur.qn_n;
        MAP.model_amp = cur.model_amp;

        MAP2.kappa = cur2.kappa;
        MAP2.sigma = cur2.sigma;
        MAP2.phi   = cur2.phi;
        MAP2.Q     = cur2.Q;
        MAP2.t0    = cur2.t0;
        MAP2.qn_n  = cur2.qn_n;
        MAP2.model_amp = cur2.model_amp;
    end

    % Thinning storage (for diagnostic plots / posterior summaries)
    if mod(k, thin_step) == 0
        ii = k / thin_step;
        t0_store(:,ii)   = cur.t0(:);
        t0_store2(:,ii)  = cur2.t0(:);
        qn_n_store(:,ii) = cur.qn_n(:);
        qn_n_store2(:,ii)= cur2.qn_n(:);
    end

    % Periodic diagnostics and saving
    if mod(k, 20000) == 0 || k == 100 || k == 1000
        % --- Posterior histograms ---
        figure(13); clf
        [~,edges] = histcounts(log10(10.^kappastore(1:k)), 50);
        histogram(10.^kappastore(1:k), 10.^edges, 'Normalization','Probability');
        set(gca,'XScale','log'); xlim([1e-10 1e-6]); xticks([1e-9 1e-8 1e-7]);
        xline(10^MAP.kappa,'r'); xlabel('Permeability (m^2)'); title('Sample 1: permeability');
        set(gca,'YTickLabel',[]);

        figure(14); clf
        histogram(sigmastore(1:k),'Normalization','Probability','NumBins',50);
        xlim([0 0.05]); xticks([0 0.02 0.04]);
        xline(MAP.sigma,'r'); xlabel('\sigma (s)'); title('Sample 1: impulse timescale');
        set(gca,'YTickLabel',[]);

        figure(17); clf
        histogram(100*phistore(1:k),'Normalization','Probability','NumBins',50);
        xlim([0 100]); xline(100*MAP.phi,'r');
        xlabel('Porosity (%)'); title('Sample 1: porosity'); set(gca,'YTickLabel',[]);

        figure(18); clf
        [~,edges] = histcounts(log10(10.^kappastore2(1:k)), 50);
        histogram(10.^kappastore2(1:k), 10.^edges, 'Normalization','Probability');
        set(gca,'XScale','log'); xlim([1e-10 1e-6]); xticks([1e-9 1e-8 1e-7]);
        xline(10^MAP2.kappa,'r'); xlabel('Permeability (m^2)'); title('Sample 2: permeability');
        set(gca,'YTickLabel',[]);

        figure(19); clf
        histogram(gasstore2(1:k)/(1.2e-6), 'Normalization','Probability','NumBins',50);
        xlim([0 200]); xlabel('Gas flow rate (ml/s)'); title('Sample 2: gas flow');
        set(gca,'YTickLabel',[]);

        figure(21); clf
        histogram(gasstore(1:k)/(1.2e-6), 'Normalization','Probability','NumBins',50);
        xlim([0 200]); xlabel('Gas flow rate (ml/s)'); title('Sample 1: gas flow');
        set(gca,'YTickLabel',[]);

        figure(22); clf
        histogram(100*phistore2(1:k),'Normalization','Probability','NumBins',50);
        xlim([0 100]); xline(100*MAP2.phi,'r');
        xlabel('Porosity (%)'); title('Sample 2: porosity'); set(gca,'YTickLabel',[]);

        % --- Model fit plots ---
        lf = length(f_h);

        figure(15); clf
        plot(f_h, obs_amp_1); hold on
        plot(f_h, MAP.model_amp(1:lf)); hold off
        set(gca,'FontSize',18); xlim([0 60]);
        xlabel('Frequency (Hz)'); ylabel('Amplitude (Pa*s)');
        title('Sample 1: observed vs MAP model');

        figure(20); clf
        plot(f_h, obs_amp_2); hold on
        plot(f_h, MAP2.model_amp(1:lf)); hold off
        set(gca,'FontSize',18); xlim([0 60]);
        xlabel('Frequency (Hz)'); ylabel('Amplitude (Pa*s)');
        title('Sample 2: observed vs MAP model');

        % --- Acceptance ratios by move type ---
        accept_by_move = zeros(1,8);
        for j = 1:8
            accept_by_move(j) = mean(acceptance(j:8:k));
        end

        % Save a compact results file (adjust as needed)
        save(sprintf('expinvresult_%d.mat', kk), ...
            'kk','tau','dt','NN','niter','window_size','idx', ...
            'kappastore','sigmastore','gasstore','phistore', ...
            'kappastore2','sigmastore2','gasstore2','phistore2', ...
            'llstore','acceptance','accept_by_move', ...
            't0_store','t0_store2','qn_n_store','qn_n_store2', ...
            'MAP','MAP2', '-v7.3');
    end

end

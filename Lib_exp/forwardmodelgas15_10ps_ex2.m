function [Posc, Pspect, A_res, A_exc] = forwardmodelgas15_10ps_ex2(m, t0, qn, sigma_k)
%FORWARDMODELGAS15_10PS_EX2 Forward model for gas-pocket pressure spectrum.
%
%   [Posc, Pspect, A_res, A_exc] = forwardmodelgas15_10ps_ex2(m, t0, qn, sigma_k)
%
% Purpose
%   Computes the steady-state pressure response of a leaky gas pocket driven
%   by a sequence of Gaussian gas-mass impulses. The model is evaluated in
%   the frequency domain and then inverse-FFT is
%   used to reconstruct a time-domain pressure series.
%
% Inputs
%   m        : parameter vector
%              m = [T, Q, R, L, kappa, D, phi]
%              T     (K)     temperature
%              Q     (kg/s)  mean gas mass flux
%              R     (m)     conduit radius
%              L     (m)     porous length
%              kappa (m^2)   Darcian permeability
%              D     (m)     gas pocket thickness
%              phi   (-)     porosity
%   t0       : 1xN_exc impulse times (s), within [0, tau]
%   qn       : 1xN_exc gas mass per impulse (kg)
%   sigma_k  : impulse duration (s), Gaussian timescale
%
% Outputs
%   Posc   : 1xN time-domain pressure (Pa), offset to oscillatory component
%   Pspect : 1xN complex pressure spectrum A_p(omega) (positive-freq stored)
%   A_res  : 1xN complex source/oscillator response function
%   A_exc  : 1xN complex excitation spectrum (sum of Gaussian impulses)
%


%% ----------------------------- Validation ------------------------------
% Basic input checks (kept lightweight for speed in MCMC).
validateattributes(m, {'numeric'}, {'vector','numel',7,'real','finite'});
validateattributes(t0, {'numeric'}, {'vector','real','finite','nonnegative'});
validateattributes(qn, {'numeric'}, {'vector','real','finite'});
validateattributes(sigma_k, {'numeric'}, {'scalar','real','finite','nonnegative'});

t0 = t0(:).';   % ensure row
qn = qn(:).';   % ensure row
assert(numel(t0) == numel(qn), 't0 and qn must have the same length.');

global tau
assert(~isempty(tau) && isscalar(tau) && tau > 0, 'Global tau must be a positive scalar.');

%% -------------------------- Physical constants --------------------------
mu_g = 1e-5;      % (Pa*s) gas viscosity
M    = 0.029;     % (kg/mol) molecular weight (air-like; update if needed)
Rg   = 8.3145;    % (J/mol/K) ideal gas constant
Pex  = 101325;    % (Pa) external/atmospheric pressure

%% --------------------------- Unpack parameters --------------------------
T     = m(1);     % (K)
Q     = m(2);     % (kg/s)
R     = m(3);     % (m)
L     = m(4);     % (m)
kappa = m(5);     % (m^2)
D     = m(6);     % (m)
phi   = m(7);     % (-)

S = pi * R^2;     % (m^2) cross-sectional area

%% -------------------------- Auxiliary parameters -------------------------
% Denominator term that appears repeatedly
den = (Pex - (Rg*T*Q^2) / (S^2 * phi^2 * M * Pex));

beta_a = S * phi * M / (Rg * T);
beta_b = mu_g * phi / (kappa * den);
beta_c = Pex * M / (Rg * T * den);
beta_d = 2 * Q / (S * phi * den);
beta_e = S * M * D / (Rg * T);

% Baseline pressure offset (P0) used to center the oscillatory pressure
P0 = Pex + (mu_g * Rg * T * Q * L) / (S * kappa * M * den);

%% --------------- Harmonic-oscillator (source) coefficients --------------
GAMMA0 = 1;
GAMMA1 = (2*(beta_a*beta_d + beta_b*beta_e)*L + beta_a*beta_b*L^2) / (2*beta_a);
GAMMA2 = (2*(beta_c*beta_e)*L + beta_a*beta_c*L^2) / (2*beta_a);

gamma0 = beta_b * L / beta_a;
gamma1 = beta_c * L / beta_a;

%% -------- Natural frequency and critical thickness (computed, optional) -
% Kept for completeness; not used in outputs.
% fn = sqrt((sqrt((GAMMA2*gamma0^2 + gamma1^2)^2 - GAMMA1^2*gamma0^2*gamma1^2) ...
%      - GAMMA2*gamma0^2) / (GAMMA2*gamma1^2)) / (2*pi);
%
% a0 = (beta_a*beta_b*L^2 + 2*beta_a*beta_d*L)^2*beta_b^2 ...
%    - 4*beta_a^2*beta_b^2*beta_c*L^2 - 4*beta_a^2*beta_c^2;
% a1 = 4*(beta_a*beta_b*L^2 + 2*beta_a*beta_d*L)*beta_b^3*L ...
%    - 8*beta_a*beta_b^2*beta_c*L;
% a2 = 4*beta_b^4*L^2;
% Dcrit = (Rg*T/(S*M))*(-a1 + sqrt(a1^2 - 4*a0*a2)) / (2*a2);

%% ------------------------ Frequency grid setup --------------------------
% Sampling rate (Hz) used for spectral grid; chosen to match sensor.
max_freq = 5000;
dt = 1 / max_freq;

time = 0:dt:tau;
N    = numel(time);

df = (1/dt) / N;  %#ok<NASGU>  % retained for consistent scaling in ifft line

% We compute only non-negative frequencies up to Nyquist:
f = 0:(1/tau):(1/(2*dt));   % Hz
Nf = numel(f);

%% -------------------- Frequency-domain response -------------------------
A_res = complex(zeros(1, Nf));
A_exc = complex(zeros(1, Nf));
A_p   = complex(zeros(1, Nf));

i = 1i;  %#ok<NASGU> keep explicit for readability

for mm = 1:Nf
    ome = 2*pi*f(mm);

    % Source/oscillator response (transfer function)
    num = (gamma0*(1i*ome)^0 + gamma1*(1i*ome)^1);
    den_tf = (GAMMA0*(1i*ome)^0 + GAMMA1*(1i*ome)^1 + GAMMA2*(1i*ome)^2);
    A_res(mm) = num / den_tf;

    % Excitation spectrum: sum of Gaussian-smoothed impulses
    %   exp(-i*omega*t0) is time shift
    %   exp(-0.5*sigma^2*omega^2) is Gaussian smoothing in time
    A_exc(mm) = sum(qn .* exp(-1i*ome*t0) .* exp(-0.5*(sigma_k^2)*ome^2));

    % Pressure spectrum
    A_p(mm) = A_res(mm) * A_exc(mm);
end

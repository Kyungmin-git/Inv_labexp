function candidatesigma = generate_sig2(currentsigma)
%GENERATE_SIG2 Propose impulse timescale sigma for Metropolis–Hastings MCMC.
%
%   candidatesigma = generate_sig2(currentsigma)
%
% Purpose
%   Random-walk proposal for the gas-impulse duration (sigma), constrained
%   to a physically/plausibly allowed interval using reflective boundaries.
%   Reflection (rather than clipping) preserves proposal symmetry near the
%   bounds, which is desirable for Metropolis–Hastings.
%
% Input
%   currentsigma   : current sigma value (s)
%
% Output
%   candidatesigma : proposed sigma value (s), guaranteed within [sigmaMin, sigmaMax]
%
% Globals (set in main script)
%   Stepsize1      : proposal standard deviation for sigma (s)
%

% -------------------------------------------------------------------------

global Stepsize1

% --- Input validation (lightweight) ---
validateattributes(currentsigma, {'numeric'}, {'scalar','real','finite'});

% --- Bounds (seconds) ---
sigmaMin = 0.005;
sigmaMax = 0.05;

% --- Random-walk proposal ---
candidatesigma = currentsigma + Stepsize1 * randn();

% --- Reflect into bounds ---
candidatesigma = reflectToBounds(candidatesigma, sigmaMin, sigmaMax);

end

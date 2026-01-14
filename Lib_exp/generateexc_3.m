function candidatet0 = generateexc_3(currentt0)
%GENERATEEXC_3 Propose new impulse-time vector t0 for MCMC (random-walk).
%
%   candidatet0 = generateexc_3(currentt0)
%
% Purpose
%   Proposes a candidate set of impulse times by perturbing ONE randomly
%   chosen element of t0 with a Gaussian step, then enforcing bounds
%   within [0, tau] using reflection at the boundaries.
%
% Inputs
%   currentt0    : 1xN (or Nx1) vector of impulse times (s), expected within [0, tau]
%
% Outputs
%   candidatet0  : same size as currentt0, with one element perturbed and
%                 all elements guaranteed to lie in [0, tau]
%

% -------------------------------------------------------------------------

global tau
assert(~isempty(tau) && isscalar(tau) && tau > 0, 'Global tau must be a positive scalar.');

% Ensure row vector (preserve original shape at the end)
origSize = size(currentt0);
t0 = currentt0(:).';  % row

N = numel(t0);
assert(N >= 1, 'currentt0 must be non-empty.');

% Choose one impulse time to perturb
idx = randi(N);

% Gaussian random-walk step (seconds)
stepStd = 0.03;            % (s) proposal standard deviation
t0(idx) = t0(idx) + stepStd * randn();

% Reflect into [0, tau]
% First reflect negative values about 0
t0(idx) = abs(t0(idx));

% Then reflect values above tau back into the interval
if t0(idx) > tau
    t0(idx) = 2*tau - t0(idx);
end

% Numerical safety: if step was huge, double-reflection could still leave it out.
% Wrap reflection until within bounds (rare, but safe).
while t0(idx) < 0 || t0(idx) > tau
    if t0(idx) < 0
        t0(idx) = -t0(idx);
    elseif t0(idx) > tau
        t0(idx) = 2*tau - t0(idx);
    end
end

% Return with original shape
candidatet0 = reshape(t0, origSize);

end

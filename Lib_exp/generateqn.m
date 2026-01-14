function [candidateqn_n, candidatelpqn_n] = generateqn(currentqn_n)
%GENERATEQN Propose impulse-weight vector qn_n and compute its log-prior.
%
%   [candidateqn_n, candidatelpqn_n] = generateqn(currentqn_n)
%
% Purpose
%   Proposes a new set of nonnegative impulse weights (dimensionless ratios)
%   used to distribute gas mass among NN impulses. The proposal perturbs ONE
%   randomly chosen element with a Gaussian step, enforces nonnegativity,
%   renormalizes so the weights sum to 1, and returns the associated
%   Gaussian log-prior that penalizes deviation from uniform weights.
%
% Input
%   currentqn_n      : 1xNN (or NNx1) current weight vector, ideally sum = 1
%
% Outputs
%   candidateqn_n    : proposed weight vector (same size as input), sum = 1
%   candidatelpqn_n  : log-prior value for candidateqn_n
%
% Prior
%   Uses a Gaussian penalty around uniform weights:
%       qn_n_mean = (1/NN) * ones(1,NN)
%       log p(q)  = - sum((q - qn_n_mean).^2) / (2*sigma_q^2)
%   with sigma_q = (0.05 * 0.2), consistent with your original code.
%


% Ensure row vector and preserve original shape at the end
origSize = size(currentqn_n);
q = currentqn_n(:).';  % row

NN = numel(q);
assert(NN >= 1, 'currentqn_n must be non-empty.');

% Choose one component to perturb
idx = randi(NN);

% Random-walk perturbation (matches original scaling: 0.05*normrnd(0,0.08))
stepStd = 0.05 * 0.08;
q(idx) = q(idx) + stepStd * randn();

% Enforce nonnegativity (keeps interpretation as weights)
q = abs(q);

% Renormalize to sum to 1 (avoid division-by-zero)
s = sum(q);
if s <= 0
    % Extremely unlikely unless everything is ~0; fallback to uniform
    q = (1/NN) * ones(1,NN);
else
    q = q / s;
end

% Compute log-prior around uniform weights
qn_n_mean = (1/NN) * ones(1,NN);
sigma_q   = (0.05 * 0.2);
candidatelpqn_n = -sum((q - qn_n_mean).^2) / (2 * sigma_q^2);

% Return with original shape
candidateqn_n = reshape(q, origSize);

end

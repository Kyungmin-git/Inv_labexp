function [candidatekappa, can_phi, can_gas] = generategas_3(currentkappa, cur_phi, cur_gas)
%GENERATEGAS_3 Propose (kappa, phi, Q) for Metropolis–Hastings MCMC.
%
%   [candidatekappa, can_phi, can_gas] = generategas_3(currentkappa, cur_phi, cur_gas)
%
% Purpose
%   Random-walk proposal for three parameters with reflective boundary handling:
%     - candidatekappa : log10(permeability [m^2])  (bounded)
%     - can_phi        : porosity (dimensionless)   (bounded)
%     - can_gas        : gas mass flux Q (kg/s)     (bounded)
%
% Bound handling
%   Uses reflection at the boundaries (instead of clipping) to maintain a
%   symmetric proposal near bounds, which is desirable for MH MCMC.
%
% Global step sizes (set in the main script)
%   Stepsize2 : step std for kappa
%   Stepsize4 : step std for phi
%   Stepsize5 : step std for gas flux
%
% Parameter bounds (edit if needed)
%   kappa in [ -11,  -6.3 ]    (log10 m^2)
%   phi   in [ 1e-3, 1.0  ]
%   gas   in [ 0,    6e-4 ]    (kg/s)  ~ 0–500 ml/s for 1.2e-6 conversion
%
% -------------------------------------------------------------------------

global Stepsize2 Stepsize4 Stepsize5

% --- Input validation (lightweight) ---
validateattributes(currentkappa, {'numeric'}, {'scalar','real','finite'});
validateattributes(cur_phi,      {'numeric'}, {'scalar','real','finite'});
validateattributes(cur_gas,      {'numeric'}, {'scalar','real','finite'});

% --- Bounds ---
kappaMin = -11;
kappaMax = -6.3;

phiMin   = 1e-3;
phiMax   = 1.0;

gasMin   = 0.0;
gasMax   = 6e-4;

% --- Random-walk proposals ---
candidatekappa = currentkappa + Stepsize2 * randn();
can_phi        = cur_phi      + Stepsize4 * randn();
can_gas        = cur_gas      + Stepsize5 * randn();

% --- Reflect into bounds (symmetric handling) ---
candidatekappa = reflectToBounds(candidatekappa, kappaMin, kappaMax);
can_phi        = reflectToBounds(can_phi,        phiMin,   phiMax);
can_gas        = reflectToBounds(can_gas,        gasMin,   gasMax);

end

% ============================ Local helper ===============================
function x = reflectToBounds(x, xmin, xmax)
%REFLECTTOBOUNDS Reflect x into [xmin, xmax].
%
% Reflection preserves proposal symmetry near boundaries.
% This version is robust even if a rare large step overshoots by more than
% one interval width.

w = xmax - xmin;
if w <= 0
    error('Invalid bounds: xmax must be > xmin.');
end

% Map to [0, 2w) then reflect second half back
y = mod(x - xmin, 2*w);   % in [0, 2w)
if y > w
    y = 2*w - y;
end
x = xmin + y;

% Numerical safety
x = min(max(x, xmin), xmax);
end

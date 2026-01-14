function [candidatekappa, can_phi, can_gas] = generategas_6(currentkappa, cur_phi, cur_gas)
%GENERATEGAS_6 Propose (kappa, phi, Q) for Metropolis–Hastings MCMC (sample #2).
%
%   [candidatekappa, can_phi, can_gas] = generategas_6(currentkappa, cur_phi, cur_gas)
%
% Purpose
%   Random-walk proposal for three parameters with reflective boundary handling:
%     - candidatekappa : log10(permeability [m^2])  (bounded)
%     - can_phi        : porosity (dimensionless)   (bounded)
%     - can_gas        : gas mass flux Q (kg/s)     (bounded)
%
% Differences vs generategas_3
%   Uses Stepsize6 (instead of Stepsize4) for porosity proposal width, to
%   allow different tuning for sample #2.
%
% Globals (set in your main script)
%   Stepsize2 : step std for kappa proposal (log10 m^2)
%   Stepsize6 : step std for phi proposal
%   Stepsize5 : step std for gas flux proposal (kg/s)
%
% Bounds (edit if needed)
%   kappa in [ -11,  -6.3 ]    (log10 m^2)
%   phi   in [ 1e-3, 1.0  ]
%   gas   in [ 0,    6e-4 ]    (kg/s)
%
% Requires
%   reflectToBounds.m in the MATLAB path (same folder is fine).
%
% -------------------------------------------------------------------------

global Stepsize2 Stepsize6 Stepsize5

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
can_phi        = cur_phi      + Stepsize6 * randn();
can_gas        = cur_gas      + Stepsize5 * randn();

% --- Reflect into bounds (symmetric handling) ---
candidatekappa = reflectToBounds(candidatekappa, kappaMin, kappaMax);
can_phi        = reflectToBounds(can_phi,        phiMin,   phiMax);
can_gas        = reflectToBounds(can_gas,        gasMin,   gasMax);

end

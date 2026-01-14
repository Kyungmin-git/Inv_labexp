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

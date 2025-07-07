function dSSF = aux_fcn_dSSF(x, U, L, k)
% Derivative of the smooth saturation function SSF
% Inputs:
%   x : variable (scalar or vector)
%   U : upper limit
%   L : lower limit
%   k : smoothing factor

    eps_mu = 1e-6;  % Avoid division by near-zero mu
    mu = max((U - L) / 2, eps_mu);  % clip mu
    alpha = (U + L) / 2;

    xi = (x - alpha) / mu;
    base = 1 + xi.^k;

    
    base(base < eps_mu) = eps_mu;

    dSSF = 1 ./ (base.^(1 + 1/k));

    % Sanity check
    dSSF(~isfinite(dSSF)) = 0;
end

function [t0] = mttf_tddb(Tc)
% Calculate mean-time-to-failure of TDDB.
% Args:
%   Tc: core temperature in Kelvin
%
% Return:
%   t0: estimated lifetime in years

% required constants
Ao = 4.0;                       % empirical value
c = 1.1949e3 / (6.022 * 1.38);  % Ea/k

t0 = Ao  * exp(c / Tc);
%fprintf("%f\n", t0);
end


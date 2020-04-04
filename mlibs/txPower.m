function [Ptx] = txPower(d, Pto)
% Calculate the average transmission power based on transmission distance.
% Args:
%   d: transmission distance in km
%   Pto: constant in transmission power
%
% Return:
%   Ptx: average transmission power in W
% required constants
alpha = 3.5;
beta = 1e-4;

Ptx = Pto + beta * (d ^ alpha);
end

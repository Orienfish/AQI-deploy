function [Ptx] = txPower(d)
% Calculate the average transmission power based on transmission distance.
% Args:
%   d: transmission distance in m
%
% Return:
%   Ptx: average transmission power in W
% required constants
alpha = 3.2;
Pto = 15.9*1e-3; % in W
beta = 5e-9;

Ptx = Pto + beta * (d ^ alpha);
end

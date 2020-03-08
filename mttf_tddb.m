function [t0] = mttf_tddb(V, Tc)
% Calculate mean-time-to-failure of TDDB.
% Args:
%   V: voltage in V
%   Tc: ambient temperature in Celsius
% 
% Return:
%   t0: estimated lifetime in h, longer with higher temperature
% required constants
a = 78.0;
b = -0.081;
X = 0.759;
Y = -66.8;
Z = ???;
k = 1.38064852e-23; % Boltzmann constant
c = ???;

Tk = Tc + 273.15; % convert from Celsius to Kelvin
t0 = c * (1/V)^(a-b*Tk) * exp((X+(Y/Tk)+Z*Tk)/(k*Tk)); 
end


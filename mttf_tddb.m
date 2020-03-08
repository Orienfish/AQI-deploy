function [t0] = mttf_tddb(V, Tc, P)
% Calculate mean-time-to-failure of TDDB.
% Args:
%   V: voltage in V
%   Tc: ambient temperature in Celsius
%   P: average power in W
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

% convert ambient temperature to core temperature
Tcore = temp_amb2core(Tc, P);
Tk = Tcore + 273.15; % convert from Celsius to Kelvin
t0 = c * (1/V)^(a-b*Tk) * exp((X+(Y/Tk)+Z*Tk)/(k*Tk)); 
end


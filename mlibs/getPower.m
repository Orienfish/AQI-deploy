function [avgPwr] = getPower(params, Tc)
% Calculate estimated power for a sensor.
%
% Args:
%   params:
%       (i) for transmission: distance in km(dtx), bandwidth in B/s (Btx), 
%                           data size in B (Ltx), constant power in W (Pto)
%       (ii) for receiving: bandwidth in B/s (Brx), data size in B (Lrx), 
%                           avg. recv. power in W (Prx)
%       (iii) for sensing: sensing power in W (Psen) and sensing time in s (tsen)
%       (iv) time for one period or data sampling period in s (T)
%       Vdd: supply voltage in V
%       f: frequency in Hz
%   not included in params:
%       Tc: core temperature in Kelvin
%
% Return:
%   avgPwr: average power consumption in W

% calculate avg. transmission power and time
Ptx = txPower(params.dtx, params.Pto);
ttx = params.Ltx / params.Btx;
%fprintf('tx power: %f tx time: %f\n', Ptx, ttx);

% calculate avg. receiving time
trx = params.Lrx / params.Brx;
%fprintf('rx power: %f rx time: %f\n', params.Prx, trx);

% calculate avg. SoC power
a = 1.5917e-11;
b = 8.61733e-7;
c = 0.07;
d = 0; % ignorable
Pd = a * params.Vdd * params.Vdd * params.f;
Ps = params.Vdd * (b * Tc * Tc * exp(c / Tc) + d);
PSoC = Pd + Ps;

% calculate total energy and power in one period
avgPwr = (Ptx * ttx + params.Prx * trx + params.Psen * params.tsen) / params.T;
avgPwr = avgPwr + PSoC;
%fprintf('tx power: %f, rx power: %f, sense power: %f, soc power: %f\n', ...
%    Ptx * ttx / params.T, params.Prx * trx / params.T, ...
%    params.Psen * params.tsen / params.T, PSoC);
end


function [avgPwr] = avgPower(dtx, Btx, Ltx, Pto, Brx, Lrx, Prx, Psen, tsen, Pslp, T)
% Calculate estimated power for a sensor.
%
% Args:
% (i) for transmission: distance in km(dtx), bandwidth in B/s (Btx), 
%                       data size in B (Ltx), constant power in W (Pto)
% (ii) for receiving: bandwidth in B/s (Brx), data size in B (Lrx), avg. recv. power in W (Prx)
% (iii) for sensing: sensing power in W (Psen) and sensing time in s (tsen)
% (iv) for sleep: avg. sleep power in W (Pslp)
% (v) time for one period or data sampling period in s (T)
%
% Return:
%   avgPwr: average power consumption in W

% calculate avg. transmission power and time
Ptx = txPower(dtx, Pto);
ttx = Ltx / Btx;
% fprintf('tx power: %f tx time: %f\n', Ptx, ttx);

% calculate avg. receiving time
trx = Lrx / Brx;
% fprintf('rx power: %f rx time: %f\n', Prx, trx);

% calculate total energy and power in one period
energy = Ptx * ttx + Prx * trx + Psen * tsen + Pslp * (T - ttx - trx - tsen);
% fprintf('tx energy: %f rx energy: %f sleep energy: %f\n', ...
%    Ptx * ttx, Prx * trx, Pslp * (T - ttx - trx));
avgPwr = energy / T;
end


function [stbPwr, stbTc] = stbPower(params, Tamb)
% Calculate the stabilized power. Iteratively fit the power and core temp.
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
%       Tamb: ambient environment temperature in Celsius.
%
% Return:
%   stbPwr: stabilized power consumption in W
%   stbTc: stabilized core temperature in Kelvin
curPwr = 0;
Tc_old = Tamb;    % initial core temperature in Celsius
eps = 1e-2;       % error bound in the iteration
iter = 0;
while 1
    curPwr = getPower(params, Tc_old + 273.15);
    Tc_new = temp_amb2core(Tamb, curPwr, Tc_old);
    if Tc_new - Tc_old < eps
        break
    end
    Tc_old = Tc_new;
    fprintf('it: %d Tc_new: %f Tc_old: %f, curPwr: %f\n', iter, ...
       Tc_new, Tc_old, curPwr);
    iter = iter + 1;
end
stbPwr = curPwr;
stbTc = Tc_old + 273.15; % convert to Kelvin
end


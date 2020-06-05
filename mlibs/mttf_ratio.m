%% Calculate mean-time-to-failure in a ratio to standard case (25 Celsius).
% Args:
%   Tc: list of core temperature in Kelvin
%
% Return:
%   MTTF: estimated MTTF ratio to standard environment of 25 Celsius

function MTTFr = mttf_ratio(Tc)
% required constants
%Ao = 4.0;                            % empirical value
c = 10 * 1.1949e3 / (6.022 * 1.38);  % Ea/k

% referece temperature (25 Celsius) in Kelvin
%Tref = 25;          % Celsius;
%Iref = 50;          % mA
%Pref = 3.3 * Iref / 1000; % reference power in W
%Tc_old = Tref;
%eps = 1e-1;       % error bound in the iteration
%iter = 0;
%while 1
%    Tc_new = amb2core(Tref, Pref, Tc_old);
%    if Tc_new - Tc_old < eps
%        break;
%    end
%    Tc_old = Tc_new;
%    iter = iter + 1;
%    fprintf('it: %d Tc_new: %f Tc_old: %f, curPwr: %f\n', iter, ...
%       Tc_new, Tc_old, Pref);
%end
%Tcref = Tc_new;
%fprintf('Tcref in Celsius: %f\n', Tcref);
Tcref = 33.010122 + 273.15; % core temperature reference in Kelvin obtained from 
                     % the code above

MTTFr = exp(c ./ Tc) / exp(c / Tcref);
%fprintf("%f\n", MTTFr);
end

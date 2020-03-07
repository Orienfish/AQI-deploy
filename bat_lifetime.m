function [t0] = bat_lifetime(y0, T, I, dt)
% Calculate battery lifetime.
% Args:
%   y0: total initial capacity in mAh
%   T: ambient temperature in Celsius
%   I: average current draw in mA
%   dt: returned lifetime resolution in h
% 
% Return:
%   t0: estimated lifetime in h, longer with higher temperature
% required constants
c = 0.56418;
Ea = 1.1949;
A = 0.96397;
R = 0.008314;

% correct initial capacity according to ambient temperature
% valid range for temperature is -5 to 40 Celsius
CF = 1.0;
if T >= -5.0 && T < 10.0
    CF = (-5.117e-7)*(T+5)^3 + 1.0076e-3*(T+5) + 9.98e-1;
elseif T >= 10.0 && T < 25.0
    CF = 2.2375e-6*(T-10)^3 + (-2.3027e-5)*(T-10)^2 + 6.6620e-4*(T-10)^1 + 1.0114;
elseif T >= 25.0 && T < 32.5
    CF = (-2.0925e-5)*(T-25)^3 + 7.7663e-5*(T-25)^2 + 1.4817e-3*(T-25)^1 + 1.0237;
elseif T >= 32.5 && T < 40
    CF = 1.7473e-5*(T-32.5)^3 + (-3.9315e-4)*(T-32.5)^2 + (-8.8444e-4)*(T-32.5)^1 + 1.0303;
end
%fprintf("%f\n", CF);
y0 = y0 * CF;
i0 = y0 * c;
j0 = y0 * (1 - c);
k = A * exp(-Ea / (R * (T + 273.15))); % note the temperature here is in Kelvin!
%fprintf("%f, %f\n", i0, k);

% iteratively compute battery lifetime
t0 = 0.0;
while i0 > 0
    t0 = t0 + dt;
    y0 = i0 + j0;
    i0 = i0 * exp(-k*dt) + (y0*k*c-I)*(1-exp(-k*dt))/k - I*c*(k*dt-1+exp(-k*dt))/k;
    j0 = j0 * exp(-k*dt) + y0*(1-c)*(1-exp(-k*dt)) - I*(1-c)*(k*dt-1+exp(-k*dt))/k;
    %fprintf("%f, %f\n", t0, i0);
end
%fprintf("%f\n", t0);
end


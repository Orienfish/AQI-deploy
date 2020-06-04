function r = bat_ratio(y, Tc, I, dt)
% Calculate battery lifetime as a ratio to the corresponding lifetime at
% 25 Celsius.
% Args:
%   y: total initial capacity in mAh
%   Tc: ambient temperature in Celsius
%   I: average current draw in mA
%   dt: lifetime resolution in h
% 
% Return:
%   r: estimated lifetime in ratio to the 25 Celsius case
% required constants
c = 0.56418;
Ea = 1.1949;
A = 0.96397;
R = 0.008314;

% correct initial capacity according to ambient temperature
% valid range for temperature is -5 to 40 Celsius
CF = 1.0;
if Tc >= -5.0 && Tc < 10.0
    CF = (-5.117e-7)*(Tc+5)^3 + 1.0076e-3*(Tc+5) + 9.98e-1;
elseif Tc >= 10.0 && Tc < 25.0
    CF = 2.2375e-6*(Tc-10)^3 + (-2.3027e-5)*(Tc-10)^2 + 6.6620e-4*(Tc-10)^1 + 1.0114;
elseif Tc >= 25.0 && Tc < 32.5
    CF = (-2.0925e-5)*(Tc-25)^3 + 7.7663e-5*(Tc-25)^2 + 1.4817e-3*(Tc-25)^1 + 1.0237;
elseif Tc >= 32.5 && Tc < 40
    CF = 1.7473e-5*(Tc-32.5)^3 + (-3.9315e-4)*(Tc-32.5)^2 + (-8.8444e-4)*(Tc-32.5)^1 + 1.0303;
end
%fprintf("%f\n", CF);

% compute the lifetime at Tc 
y0 = y * CF;
i0 = y0 * c;
j0 = y0 * (1 - c);
k = A * exp(-Ea / (R * (Tc + 273.15))); % note the temperature here is in Kelvin!
%fprintf("i0: %f, k: %f\n", i0, k);

% compute the lifetime at reference temperature of 25 Celsius
% and reference current draw of 50mA
Tref = 25;          % Celsius
CFref = 1.0237;
Iref = 50;          % mA
% compute the lifetime at Tc 
yref = y0 * CFref;
iref = yref * c;
jref = yref * (1 - c);
kref = A * exp(-Ea / (R * (Tref + 273.15))); % note the temperature here is in Kelvin!
%fprintf("iref: %f, kref: %f\n", iref, kref);

% iteratively compute battery lifetime
t0 = 0.0;
tref = 0.0;
while i0 > 0 || iref > 0
    if i0 > 0
        t0 = t0 + dt;
        y0 = i0 + j0;
        i0 = i0 * exp(-k*dt) + (y0*k*c-I)*(1-exp(-k*dt))/k - I*c*(k*dt-1+exp(-k*dt))/k;
        j0 = j0 * exp(-k*dt) + y0*(1-c)*(1-exp(-k*dt)) - I*(1-c)*(k*dt-1+exp(-k*dt))/k;
        %fprintf("%f, %f\n", t0, i0);
    end
    if iref > 0
        tref = tref + dt;
        yref = iref + jref;
        iref = iref * exp(-kref*dt) + (yref*kref*c-Iref)*(1-exp(-kref*dt))/kref ...
            - Iref*c*(kref*dt-1+exp(-kref*dt))/kref;
        jref = jref * exp(-kref*dt) + yref*(1-c)*(1-exp(-kref*dt)) - ...
            Iref*(1-c)*(kref*dt-1+exp(-kref*dt))/kref;
        %fprintf("%f, %f\n", tref, iref);
    end
end
%fprintf("t0: %f\n", t0);
%fprintf("tref: %f\n", tref);

% compute the ratio
r = t0 / tref;
%fprintf("ratio: %f\n", r);
end


function [Tcore] = temp_amb2core(Tamb, pwr, Tcore_last)
% Iteratively calculate the ambient temperature for next time stamp.
% Args:
%   Tamb: ambient temperature in Celsius
%   pwr: average power in W
%   Tcore_last: the core temperature in last time stamp in Kelvin? Celsius?
%
% Return:
%   Tcore: core processor temperature in Celsius
% required constants
m_A = 1.41558165;
m_B = 0.98;
m_C = 4.43815458;
m_D = 1.7971150984226512;

Tcore = (m_A * (1-m_B) * Tamb + m_B * Tcore_last + m_C * (1-m_B) * pwr + m_D * (1-m_B));
fprintf("%f\n", Tcore);
end


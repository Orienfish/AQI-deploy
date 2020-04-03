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
m_A = 0.0283;
m_B = 0.98;
m_C = 0.0888;
m_D = 0.0359;

Tcore = m_A * Tamb + m_B * Tcore_last + m_C * pwr + m_D;
% fprintf("%f\n", Tcore);
end


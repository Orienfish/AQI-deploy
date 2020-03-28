function [Tcore] = temp_amb2core(Tamb, pwr)
% Convert ambient temperature to core temperature.
% Args:
%   Tamb: ambient temperature in Celsius
%   pwr: average power in W
%
% Return:
%   Tcore: core processor temperature in Celsius
% required constants
m_A = 1.41558165;
m_B = 4.43815458;
m_C = 1.7971150984226512;

Tcore = (m_A * Tamb + m_B * pwr + m_C);
%fprintf("%f\n", Tcore);
end


function [y] = Jinc(x)
%Jinc Summary of this function goes here
%   Detailed explanation goes here
% Jinc(x) = J1(2*pi*x)/x

% locate non-zero elements of x:
mask = (x ~= 0);
% initialize output with pi:
y = pi*ones(size(x));
% compute output values for all other x:
y(mask) = besselj(1,2*pi*x(mask))./(x(mask));
end
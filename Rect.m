function y = Rect(v)
%RECT Summary of this function goes here
%   Detailed explanation goes here
%y = H(v+.5) - H(v-.5);
y = abs(v/2) <= 1/2;
end


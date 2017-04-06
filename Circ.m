function [z] = Circ(x,y)
%Circ Summary of this function goes here
%   Detailed explanation goes here
z = 0+(sqrt((x/2).^2 + (y/2).^2) <= 1/2);

end


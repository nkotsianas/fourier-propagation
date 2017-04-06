function y = H(v)
%H Summary of this function goes here
%   Detailed explanation goes here

y = zeros(1,length(v));
for k = 1:length(v)
    if v(k) == 0
        y(k) = .5;
    elseif v(k) > 0
        y(k) = 1;
    else
        y(k) = 0;
    end
end

end


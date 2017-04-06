function s = Step(fHandle,x,N)
%STEP Summary of this function goes here
%   Detailed explanation goes here
s = fHandle(floor(N*x)/N);
end
% Convolution of two Gaussian functions:

wa = .3;    % Gaussian 1 width [exp(-pi) radius] (m)
wb = .2;    % Gaussian 2 width [exp(-pi) radius] (m)
L = 2;      % side length (m)
M = 200;    % number of samples
dx = L/M;   % sample interval (m)

x = -L/2:dx:L/2-dx;             % x coordinates
fa = exp(-pi*(x.^2)/(wa^2));    % Gaussian a
fb = exp(-pi*(x.^2)/(wb^2));    % Gaussian b

figure(1)
plot(x,fa,x,fb,'--'); title('functions'); xlabel('x (m)');

Fa = dx*ifftshift(fft(fftshift(fa)));           % fft of fa
Fb = dx*ifftshift(fft(fftshift(fb)));           % fft of fb
FaFb = Fa.*Fb;                                  % product of the functions
fafb = 1/dx*fftshift(ifft(ifftshift(FaFb)));    % inverse transform, scale, and center

figure(2)
plot(x,fafb); title('convolution'); xlabel('x (m)');
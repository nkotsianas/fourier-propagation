% 2D FFT example:

wx = .015;
wy = .015;
Lx = .2;
Ly = .2;
M = 200;
N = 200;
dx = Lx/M;
dy = Ly/N;
x = -Lx/2:dx:Lx/2-dx;
y = -Ly/2:dy:Ly/2-dy;
[X,Y] = meshgrid(x,y);

g = Circ(X/(2*wx),Y/(2*wy));

figure(1)
imagesc(x,y,g);     % image display
colormap('gray');   % linear gray display map
axis square;        % sqare figure
axis xy;            % y positive up
xlabel('x (m)'); ylabel('y (m)');

figure(2)
plot(x,g(M/2+1,:));     % slice profile
xlabel('x (m)'); ylabel('g(M/2+1,x)');

gs = fftshift(g);
Gs = dx*dy*fft2(gs);
G = ifftshift(Gs);

Lfx = 1/dx;
Lfy = 1/dy;
dfx = 1/Lx;
dfy = 1/Ly;
fx = -Lfx/2:dfx:Lfx/2-dfx;
fy = -Lfy/2:dfy:Lfy/2-dfy;


figure(8)
surf(x,y,g)
camlight left; lighting phong
colormap('gray')
shading interp

figure(3)
surf(fx,fy,abs(G))
camlight left; lighting phong
colormap('gray')
shading interp
xlabel('fx (cyc/m)'); ylabel('fy (cyc/m)');

figure(4)
plot(fx,abs(G(M/2+1,:)));
title('magnitude of G'); xlabel('fx (cyc/m)');












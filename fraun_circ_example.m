% Fraunhofer irradiance plot:

Lx = .2;
Ly = .2;
M = 250;
N = 250;
dx = Lx/M;
dy = Ly/N;
x = -Lx/2:dx:Lx/2-dx;
y = -Ly/2:dy:Ly/2-dy;
[X,Y] = meshgrid(x,y);

R = 1e-3;
lambda = 633e-9;
z = 50;
k = 2*pi/lambda;
lz = lambda*z;

% irradiance:
U2 = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2))*R^2.*Jinc(R/lz*sqrt(X.^2 + Y.^2));
I2 = abs(U2).^2;
% I2 = (R^2/lz)^2.*(Jinc(R/lz*sqrt(X.^2 + Y.^2))).^2;

% figure(1)
% imagesc(x,y,U1);
% colormap('gray');
% axis square;
% axis xy;
% 
% return;

figure(1)   % irradiance image
imagesc(x,y,I2);
colormap('gray');
axis square;
axis xy;

figure(2)
imagesc(x,y,nthroot(I2,3));
colormap('gray');
axis square;
axis xy;

figure(3)   % x-axis profile
plot(x,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');

figure(4)
surf(X,Y,I2);
camlight left; lighting phong
colormap('jet')
shading interp



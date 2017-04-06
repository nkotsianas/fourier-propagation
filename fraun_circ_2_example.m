% Fraunhofer irradiance plot of circular aperture:

R = 1e-3;
lambda = 633e-9;
z = 50;
k = 2*pi/lambda;
lz = lambda*z;


Lx = .2;
Ly = .2;
Mx = 256;
Ny = 256;
dx = Lx/Mx;
dy = Ly/Ny;
x = -Lx/2:dx:Lx/2-dx;
y = -Ly/2:dy:Ly/2-dy;
[X,Y] = meshgrid(x,y);

Lfx = 1/dx;
Lfy = 1/dy;
dfx = 1/Lx;
dfy = 1/Ly;
fx = -Lfx/2:dfx:Lfx/2-dfx;
fy = -Lfy/2:dfy:Lfy/2-dfy;
[fX,fY] = meshgrid(fx,fy);



Lxi = lz/dx;
Leta = lz/dy;
Mxi = Mx;
Neta = Ny;
dxi = Lxi/Mxi;
deta = Leta/Neta;
xi = -Lxi/2:dxi:Lxi/2-dxi;
eta = -Leta/2:deta:Leta/2-deta;
[XI,ETA] = meshgrid(xi,eta);


Lfxi = 1/dxi;
Lfeta = 1/deta;
dfxi = 1/Lxi;
dfeta = 1/Leta;
fxi = -Lfxi/2:dfxi:Lfxi/2-dfxi;
feta = -Lfeta/2:dfeta:Lfeta/2-dfeta;
[fXI,fETA] = meshgrid(fxi,feta);



% development of irradiance plot
U1 = Circ(XI/(2*R),ETA/(2*R));
FU1 = dxi*deta*ifftshift(fft2(fftshift(U1)));
FU1A = R^2*Jinc(R*sqrt(fXI.^2 + fETA.^2));

FU1known = R^2*Jinc(R/lz*sqrt(X.^2 + Y.^2));

figure(1)
surf(fXI,fETA,abs(FU1));
camlight left; lighting phong
colormap('jet')
shading interp
xlabel('fXI'); ylabel('fETA');

figure(2)
surf(fXI,fETA,abs(FU1A));
camlight left; lighting phong
colormap('jet')
shading interp
xlabel('fXI'); ylabel('fETA');

figure(3)
surf(fX,fY,abs(FU1known));
camlight left; lighting phong
colormap('jet')
shading interp
xlabel('fX'); ylabel('fY');
% return;



U2 = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2)).*FU1;
U2A = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2)).*FU1A;
U2known = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2)).*FU1known;

I2 = abs(U2).^2;
I2A = abs(U2A).^2;
I2known = abs(U2known).^2;

figure(4)   % irradiance plot
imagesc(x,y,nthroot(I2,3));
colormap('jet');
axis square;
axis xy;
figure(5)
imagesc(x,y,nthroot(I2A,3));
colormap('jet');
axis square;
axis xy;
figure(6)
imagesc(x,y,nthroot(I2known,3));
colormap('jet');
axis square;
axis xy;


figure(7)   % x-axis profile
plot(x,I2(Ny/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
figure(8)
plot(x,I2A(Ny/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
figure(9)
plot(x,I2known(Ny/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');







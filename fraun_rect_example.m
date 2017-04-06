% Fraunhofer irradiance plot of rectangular aperture:

wxi = .1e-3;
weta = .05e-3;
lambda = 633e-9;
z = 5;
k = 2*pi/lambda;
lz = lambda*z;


Lx = .3;
Ly = .6;
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
U1 = Rect(XI/(2*wxi)).*Rect(ETA/(2*weta));
FU1 = dxi*deta*ifftshift(fft2(fftshift(U1)));
FU1A = 4*wxi*weta*sinc(2*wxi*fXI).*sinc(2*weta*fETA);


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

% return;



U2 = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2)).*FU1;
U2A = exp(1i*k*z)/(1i*lz)*exp(1i*k/(2*z)*(X.^2 + Y.^2)).*FU1A;

I2 = abs(U2).^2;
I2A = abs(U2A).^2;

figure(3)   % irradiance plot
imagesc(x,y,nthroot(I2,3));
colormap('jet');
axis square;
axis xy;
figure(4)
imagesc(x,y,nthroot(I2A,3));
colormap('jet');
axis square;
axis xy;


figure(5)   % x-axis profile
plot(x,I2(Ny/2+1,:),x,I2A(Ny/2+1,:));
legend('I2','I2A');
axis([x(1) x(length(x)) 0 max(I2(Ny/2+1,Mx/2+1),I2A(Ny/2+1,Mx/2+1))]);
% axis([x(1) x(length(x)) 0 1e-5]);
xlabel('x (m)'); ylabel('Irradiance');
% figure(6)
% plot(x,I2A(Ny/2+1,:));
% axis([x(1) x(length(x)) 0 max(I2A(Ny/2+1,Mx/2+1))]);
% % axis([x(1) x(length(x)) 0 1e-5]);
% xlabel('x (m)'); ylabel('Irradiance');



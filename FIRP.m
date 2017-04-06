function [U_2] = FIRP(U_1,Lx_2,Ly_2,lambda,z)
%Fresnel Impulse Response Propagator Summary of this function goes here
%   Detailed explanation goes here
k = 2*pi/lambda;
lz = lambda*z;

[Mx_1,Ny_1] = size(U_1);

% Lx_2 = .1;
% Ly_2 = .1;
Mx_2 = Mx_1;
Ny_2 = Ny_1;
dx_2 = Lx_2/Mx_2;
dy_2 = Ly_2/Ny_2;
x_2 = -Lx_2/2:dx_2:Lx_2/2-dx_2;
y_2 = -Ly_2/2:dy_2:Ly_2/2-dy_2;
[X_2,Y_2] = meshgrid(x_2,y_2);

% Lfx = 1/dx;
% Lfy = 1/dy;
% dfx = 1/Lx;
% dfy = 1/Ly;
% fx = -Lfx/2:dfx:Lfx/2-dfx;
% fy = -Lfy/2:dfy:Lfy/2-dfy;
% [fX,fY] = meshgrid(fx,fy);

Lx_1 = lz/dx_2;
Ly_1 = lz/dy_2;
% Mx_1 = Mx;
% Ny_1 = Ny;
dx_1 = Lx_1/Mx_1;
dy_1 = Ly_1/Ny_1;
% xi = -Lxi/2:dxi:Lxi/2-dxi;
% eta = -Leta/2:deta:Leta/2-deta;
% [XI,ETA] = meshgrid(xi,eta);

Lfx_1 = 1/dx_1;
Lfy_1 = 1/dy_1;
% dfxi = 1/Lxi;
% dfeta = 1/Leta;
% fxi = -Lfxi/2:dfxi:Lfxi/2-dfxi;
% feta = -Lfeta/2:dfeta:Lfeta/2-dfeta;
% [fXI,fETA] = meshgrid(fxi,feta);


h = exp(1j*k*z)/(1j*lz).*exp(1j*k/(2*z)*(X_2.^2 + Y_2.^2));
H = dx_2*dy_2*ifftshift(fft2(fftshift(h)));
FU_1 = dx_1*dy_1*ifftshift(fft2(fftshift(U_1)));
FU_2 = FU_1.*H;
U_2 = Lfx_1*Lfy_1*fftshift(ifft2(ifftshift(FU_2)));



end


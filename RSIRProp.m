function [U_2] = RSIRProp(U_1,Lx_1,Ly_1,lambda,z)
% Fraunhofer Propagator

% "_1" indicates source plane variables;
% "_2" indicates observation plane variables;

k = 2*pi/lambda;
lz = lambda*z;

[Mx_1,Ny_1] = size(U_1);
dx_1 = Lx_1/Mx_1;
dy_1 = Ly_1/Ny_1;
% x_1 = -Lx_1/2:dx_1:Lx_1/2-dx_1;
% y_1 = -Ly_1/2:dy_1:Ly_1/2-dy_1;
% [X_1,Y_1] = meshgrid(x_1,y_1);


Mx_2 = Mx_1;
Ny_2 = Ny_1;
Lx_2 = lz*Mx_2/Lx_1;
Ly_2 = lz*Ny_2/Ly_1;
dx_2 = Lx_2/Mx_2;
dy_2 = Ly_2/Ny_2;
x_2 = -Lx_2/2:dx_2:Lx_2/2-dx_2;
y_2 = -Ly_2/2:dy_2:Ly_2/2-dy_2;
[X_2,Y_2] = meshgrid(x_2,y_2);


% % Lx_2 = .1; % given as arg
% % Ly_2 = .1; % given as arg
% Mx_2 = Mx_1;
% Ny_2 = Ny_1;
% dx_2 = Lx_2/Mx_2;
% dy_2 = Ly_2/Ny_2;
% % x = -Lx/2:dx:Lx/2-dx;
% % y = -Ly/2:dy:Ly/2-dy;
% % [X,Y] = meshgrid(x,y);
% 
% % Lfx_2 = 1/dx_2;
% % Lfy_2 = 1/dy_2;
% % dfx_2 = 1/Lx_2;
% % dfy_2 = 1/Ly_2;
% % fx_2 = -Lfx_2/2:dfx_2:Lfx_2/2-dfx_2;
% % fy_2 = -Lfy_2/2:dfy_2:Lfy_2/2-dfy_2;
% % [fX_2,fY_2] = meshgrid(fx_2,fy_2);
% 
% Lx_1 = lz/dx_2;
% Ly_1 = lz/dy_2;
% % Mx_1 = Mx_2;
% % Ny_1 = Ny_2;
% dx_1 = Lx_1/Mx_1;
% dy_1 = Ly_1/Ny_1;
% % xi = -Lxi/2:dxi:Lxi/2-dxi;
% % eta = -Leta/2:deta:Leta/2-deta;
% % [XI,ETA] = meshgrid(xi,eta);
% 
Lfx_1 = 1/dx_1;
Lfy_1 = 1/dy_1;
% % dfxi = 1/Lxi;
% % dfeta = 1/Leta;
% % fxi = -Lfxi/2:dfxi:Lfxi/2-dfxi;
% % feta = -Lfeta/2:dfeta:Lfeta/2-dfeta;
% % [fXI,fETA] = meshgrid(fxi,feta);

r2 = z^2 + X_2.^2 + Y_2.^2;
h = (z/(1j*lambda))*exp(1j*k*sqrt(r2))./r2;
Fh_2 = dx_2*dy_2*ifftshift(fft2(fftshift(h)));
FU_1 = dx_1*dy_1*ifftshift(fft2(fftshift(U_1)));
FU_2 = FU_1.*Fh_2;
U_2 = Lfx_1*Lfy_1*fftshift(ifft2(ifftshift(FU_2)));

end


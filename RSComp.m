
% Direct computation of the Rayleigh-Sommerfeld diffraction integral

% all values in **microns** (for better computational accuracy)
% e3 = mm
% e4 = cm
% e6 = m
um = 1;
mm = 1e3;
cm = 1e4;

tic;
z = 30*cm;
lambda = .5*um;
k = 2*pi/lambda;

hwx_1 = 5*cm;
hwy_1 = hwx_1;
xf_1 = hwx_1;

airy_x = 2*1.22*lambda*z/(2*hwx_1);
disp(['Airy diameter in x: ', num2str(airy_x), 'um']);
xf_2 = 5*airy_x/2;
% return;

Nx_1a = 256;
Nx_2a = 64;


% xf_1 = L_1/2;
x0_1 = -xf_1;
yf_1 = xf_1;
y0_1 = -yf_1;
% Nx_1a = 100;
Ny_1a = Nx_1a;
dx_1 = (xf_1-x0_1)/(Nx_1a+1);
dy_1 = (yf_1-y0_1)/(Ny_1a+1);

x_1 = x0_1 + ((0:Nx_1a) + .5)*dx_1;
y_1 = y0_1 + ((0:Ny_1a) + .5)*dy_1;
[X_1,Y_1] = ndgrid(x_1,y_1);

Nx_1 = length(x_1);
Ny_1 = length(y_1);

% xf_2 = L_2/2;
x0_2 = -xf_2;
yf_2 = xf_2;
y0_2 = -yf_2;
% Nx_2a = 100;
Ny_2a = Nx_2a;
dx_2 = (xf_2-x0_2)/Nx_2a;
dy_2 = (yf_2-y0_2)/Ny_2a;

x_2 = x0_2 + (0:Nx_2a)*dx_2;
y_2 = y0_2 + (0:Ny_2a)*dy_2;
[X_2,Y_2] = ndgrid(x_2,y_2);

Nx_2 = length(x_2);
Ny_2 = length(y_2);


dn = .5;
f = z;
R = f*dn;

P_1 = Circ(X_1/hwx_1,Y_1/hwy_1);
% P_1 = Rect(X_1/hwx_1).*Rect(Y_1/hwy_1);

s_1 = R - sqrt(R^2 - X_1.^2 - Y_1.^2); % lens counter-thickness

phase_1 = -k*dn*s_1; % (phase_1 ignores constant 'k*n*max_thickness' term. (sorry jesus))
modphase_1 = -2*pi*Frac(dn*s_1/lambda);
tA_1 = P_1.*exp(1j*phase_1);
modtA_1 = P_1.*exp(1j*modphase_1);
% U_1 = tA_1;
U_1 = modtA_1;


% figure(1)
% imagesc(x_1,y_1,P_1);
% axis square; axis xy;
% colormap('gray');
% xlabel('x_1 (\mu{}m)'); ylabel('y_1 (\mu{}m)'); title('P_1');
% % return;
% 
% figure(2)
% plot(x_1,s_1(Ny_1a/2+1,:),x_1,s_1(Ny_1a/2+1,:).*P_1(Ny_1a/2+1,:));
% xlabel('x_1 (\mu{}m)'); title('s_1 (\mu{}m)');

% phase_1
% % figure(3)
% % imagesc(x_1,y_1,phase_1);
% % axis square; axis xy;
% % colormap('gray');
% % xlabel('x_1 (\mu{}m)'); ylabel('y_1 (\mu{}m)'); title('phase_1'); 
figure(4)
plot(x_1,phase_1(Ny_1a/2+1,:),'-o');
xlabel('x_1 (\mu{}m)'); title('phase_1: x profile');

% modphase_1
% % figure(5)
% % imagesc(x_1,y_1,modphase_1);
% % axis square; axis xy;
% % colormap('gray');
% % xlabel('x_1 (\mu{}m)'); ylabel('y_1 (\mu{}m)'); title('modphase_1'); 
figure(6)
plot(x_1,modphase_1(Ny_1a/2+1,:),'-o');
xlabel('x_1 (\mu{}m)'); title('modphase_1: x profile');

% tA_1
figure(7)
plot(x_1,real(tA_1(Ny_1a/2+1,:)),x_1,imag(tA_1(Ny_1a/2+1,:)));
legend('real','imag');
xlabel('x_1 (\mu{}m)'); title('tA_1: x profile');

% modtA_1
figure(8)
plot(x_1,real(modtA_1(Ny_1a/2+1,:)),x_1,imag(modtA_1(Ny_1a/2+1,:)));
legend('real','imag');
xlabel('x_1 (\mu{}m)'); title('modtA_1: x profile');


% % U_1
% figure(7)
% imagesc(x_1,y_1,real(U_1));
% axis square; axis xy;
% colormap('gray');
% xlabel('x_1 (\mu{}m)'); ylabel('y_1 (\mu{}m)'); title('re(U_1)'); 
% 
% figure(8)
% imagesc(x_1,y_1,imag(U_1));
% axis square; axis xy;
% colormap('gray');
% xlabel('x_1 (\mu{}m)'); ylabel('y_1 (\mu{}m)'); title('im(U_1)'); 
% figure(9)
% plot(x_1,real(U_1(Ny_1a/2+1,:)),x_1,imag(U_1(Ny_1a/2+1,:)));
% legend('real','imag');
% xlabel('x_1 (\mu{}m)'); title('U_1: x profile');

% return;


U_2 = zeros(Nx_2,Ny_2);
for kx_2 = 1:Nx_2
    for ky_2 = 1:Ny_2
        r2_12 = z^2 + (x_2(kx_2) - X_1).^2 + (y_2(ky_2) - Y_1).^2;
        r_12 = sqrt(r2_12);
        ex = 2*pi*Frac(r_12/lambda);
        exc = exp(1i*ex);
        fint = U_1.*exc./r2_12;
        U_2(kx_2,ky_2) = sum(sum(fint))*dx_1*dy_1;
    end
end
toc;

I_2 = abs(U_2).^2;

vol_I_2 = sum(sum(I_2))*dx_2*dy_2;
I_1 = abs(U_1).^2;
vol_I_1 = sum(sum(I_1))*dx_1*dy_1;

disp({'vol_I_1: ', vol_I_1;...
    'vol_I_2: ', vol_I_2;...
    'vol_I_2/vol_I_1: ', vol_I_2/vol_I_1});

figure(10)
imagesc(x_2,y_2,nthroot(I_2,3));
axis square; axis xy;
colormap('gray');
xlabel('x_2 (\mu{}m)'); ylabel('y_2 (\mu{}m)'); title('I_2'); 

figure(11)
surf(X_2,Y_2,I_2);
camlight left; lighting phong;
colormap('gray');
shading interp;
xlabel('x_2 (\mu{}m)'); ylabel('y_2 (\mu{}m)'); title('I_2'); 

figure(12)
plot(x_2,I_2(Ny_2a/2+1,:));
axis([min(x_2) max(x_2) 0 max(I_2(Ny_2a/2+1,:))]);
xlabel('x_2 (\mu{}m)'); title('I_2'); 







% Circular converging PARABOLIC Lens propagation

% "_1" indicates source plane variables;
% "_2" indicates observation plane variables;

wx_1 = 1;
wy_1 = wx_1;
lambda = 500e-9;
z = 5;

k = 2*pi/lambda;
lz = lambda*z;
% For irradiance, only sampling condition is the usual source plane...
% sampling condition: B1 <= 1/(2*dx_1);

% Critical sampling: dx_1 = lz/Lx_1;


Mx_1 = 2400000;
Ny_1 = Mx_1;
Lx_1 = 5;
Ly_1 = Lx_1;
dx_1 = Lx_1/Mx_1;
dy_1 = Ly_1/Ny_1;
x_1 = -Lx_1/2:dx_1:Lx_1/2-dx_1;
y_1 = -Ly_1/2:dy_1:Ly_1/2-dy_1;
[X_1,Y_1] = meshgrid(x_1,y_1);

Lfx_1 = 1/dx_1;
Lfy_1 = 1/dy_1;
dfx_1 = 1/Lx_1;
dfy_1 = 1/Ly_1;
fx_1 = -Lfx_1/2:dfx_1:Lfx_1/2-dfx_1;
fy_1 = -Lfy_1/2:dfy_1:Lfy_1/2-dfy_1;
[fX_1,fY_1] = meshgrid(fx_1,fy_1);

Mx_2 = Mx_1;
Ny_2 = Ny_1;
Lx_2 = lz*Mx_2/Lx_1;
Ly_2 = lz*Ny_2/Ly_1;
dx_2 = Lx_2/Mx_2;
dy_2 = Ly_2/Ny_2;
x_2 = -Lx_2/2:dx_2:Lx_2/2-dx_2;
y_2 = -Ly_2/2:dy_2:Ly_2/2-dy_2;
[X_2,Y_2] = meshgrid(x_2,y_2);

Lfx_2 = 1/dx_2;
Lfy_2 = 1/dy_2;
dfx_2 = 1/Lx_2;
dfy_2 = 1/Ly_2;
fx_2 = -Lfx_2/2:dfx_2:Lfx_2/2-dfx_2;
fy_2 = -Lfy_2/2:dfy_2:Lfy_2/2-dfy_2;
[fX_2,fY_2] = meshgrid(fx_2,fy_2);




% development of irradiance plot
f = z;
% fm = f*dn/2;
AiryZero = 3.83170597020751*f/(k*wx_1);
% display some sampling parameters
disp(' wx_1 // z // f // AiryZero (um) // Mx_1 // Lx_1');
disp([' ',num2str(wx_1),' // ',num2str(z),' // ',...
    num2str(f),' // ',num2str(AiryZero*1e6),' // ',...
    num2str(Mx_1),' // ',num2str(Lx_1)]);
disp('  >1 // TF<1<IR // Fraun<<Fres<1');
disp(['x: ', num2str(lambda*f*Mx_1/(2*wx_1*Lx_1)),' // ',...
    num2str(Lx_1^2/(Mx_1*lz)),' // ',...
    num2str(wx_1^2/lz)]);
disp(['y: ', num2str(lambda*f*Ny_1/(2*wy_1*Ly_1)),' // ',...
    num2str(Ly_1^2/(Ny_1*lz)),' // ',...
    num2str(wy_1^2/lz)]);
% return;
% P_1 = Rect(X_1/(2*wx_1)).*Rect(Y_1/(2*wy_1));
P_1 = Circ(X_1/(2*wx_1),Y_1/(2*wy_1));
% s_1 = (X_1.^2)/(4*f); % lens counter-thickness
% phase_1 = -k*dn*(X_1.^2+Y_1.^2)/(4*fm);
phase_1 = -k/(2*f)*(X_1.^2+Y_1.^2); % 'surely youre joking' phase
% (phase_1 ignores constant 'k*n*max_thickness' term. (sorry jesus))
tA_1 = P_1.*exp(1j*phase_1);
U_1 = tA_1;




% figure(2)
% surf(X_1,Y_1,phase_1);
% camlight left; lighting phong
% colormap('gray')
% shading interp
% return;
% 
figure(1)
imagesc(x_1,y_1,P_1.*phase_1);
colormap('gray');
axis square; axis xy;
xlabel('x_1 (m)'); ylabel('y_1 (m)'); title('phase_1');
figure(2)   % x-axis phase profile
plot(x_1,P_1(Ny_1/2+1,:).*phase_1(Ny_1/2+1,:));
xlabel('x_1 (m)'); title('phase_1');
% return;
% 
% 
% 
% 
% FU_1 = dx_1*dy_1*ifftshift(fft2(fftshift(U_1)));
% figure(1)
% imagesc(fx_1,fy_1,real(FU_1));
% colormap('gray');
% axis square; axis xy;
% xlabel('fx_1 (m)'); ylabel('fy_1 (m)'); title('Re[FU_1]');
% figure(2)
% imagesc(fx_1,fy_1,imag(FU_1));
% colormap('gray');
% axis square; axis xy;
% xlabel('fx_1 (m)'); ylabel('fy_1 (m)'); title('Im[FU_1]');
% figure(3)   % x-axis FU_1 profile
% plot(fx_1,real(FU_1(Ny_1/2+1,:)),fx_1,imag(FU_1(Ny_1/2+1,:)));
% xlabel('fx_1 (m)'); title('FU_1: x');
% legend('real','imag');
% figure(4)   % y-axis FU_1 profile
% plot(fy_1,real(FU_1(:,Mx_1/2+1)),fy_1,imag(FU_1(:,Mx_1/2+1)));
% xlabel('fy_1 (m)'); title('FU_1: y');
% legend('real','imag');
% 
% return;
% 


% figure(1)
% imagesc(x_1,y_1,real(U_1));
% colormap('gray');
% axis square; axis xy;
% xlabel('x_1 (m)'); ylabel('y_1 (m)'); title('Re[U_1]');
% figure(2)
% imagesc(x_1,y_1,imag(U_1));
% colormap('gray');
% axis square; axis xy;
% xlabel('x_1 (m)'); ylabel('y_1 (m)'); title('Im[U_1]');
% figure(3)   % x-axis U_1 profile
% plot(x_1,real(U_1(Ny_1/2+1,:)),x_1,imag(U_1(Ny_1/2+1,:)));
% xlabel('x_1 (m)'); title('U_1');
% legend('real','imag');
% figure(4)   % x-axis imag(U_1) profile
% plot(x_1,imag(U_1(Ny_1/2+1,:)));
% xlabel('x_1 (m)'); title('Im[U_1]');
% % return;
% 




U_2 = FraunProp(U_1,Lx_1,Ly_1,lambda,z);
UTF_2 = FTFP(U_1,Lx_2,Ly_2,lambda,z);
UIR_2 = FIRP(U_1,Lx_2,Ly_2,lambda,z);
URSIR_2 = RSIRProp(U_1,Lx_1,Ly_1,lambda,z);
URSTF_2 = RSTFProp(U_1,Lx_1,Ly_1,lambda,z);


IRSIR_2 = abs(URSIR_2).^2;
IRSTF_2 = abs(URSTF_2).^2;
I_2 = abs(U_2).^2;
ITF_2 = abs(UTF_2).^2;
IIR_2 = abs(UIR_2).^2;


Ivol_ref = Volume(abs(U_1).^2,Lx_1,Ly_1);
Ivol = Volume(I_2,Lx_2,Ly_2)/Ivol_ref;
ITFvol = Volume(ITF_2,Lx_2,Ly_2)/Ivol_ref;
IIRvol = Volume(IIR_2,Lx_2,Ly_2)/Ivol_ref;
IRSTFvol = Volume(IRSTF_2,Lx_2,Ly_2)/Ivol_ref;
IRSIRvol = Volume(IRSIR_2,Lx_2,Ly_2)/Ivol_ref;


figure(3)
imagesc(x_2,y_2,nthroot(I_2,3));
colormap('gray');
axis square; axis xy;
xlabel('x_2 (m)'); ylabel('y_2 (m)'); title(['I_2 / I_1 = ',num2str(Ivol)]);
figure(4)
plot(x_2,I_2(Ny_2/2+1,:));
xlabel('x_2 (m)'); title('I_2');

figure(5)
imagesc(x_2,y_2,nthroot(ITF_2,3));
colormap('gray');
axis square; axis xy;
xlabel('x_2 (m)'); ylabel('y_2 (m)'); title(['ITF_2 / I_1 = ',num2str(ITFvol)]);
figure(6)
plot(x_2,ITF_2(Ny_2/2+1,:));
xlabel('x_2 (m)'); title('ITF_2');

figure(7)
imagesc(x_2,y_2,nthroot(IIR_2,3));
colormap('gray');
axis square; axis xy;
xlabel('x_2 (m)'); ylabel('y_2 (m)'); title(['IIR_2 / I_1 = ',num2str(IIRvol)]);
figure(8)
plot(x_2,IIR_2(Ny_2/2+1,:));
xlabel('x_2 (m)'); title('IIR_2');





figure(9)
imagesc(x_2,y_2,nthroot(IRSTF_2,3));
colormap('gray');
axis square; axis xy;
xlabel('x_2 (m)'); ylabel('y_2 (m)'); title(['IRSTF_2 / I_1 = ',num2str(IRSTFvol)]);
figure(10)
plot(x_2,IRSTF_2(Ny_2/2+1,:));
xlabel('x_2 (m)'); title('IRSTF_2');

figure(11)
imagesc(x_2,y_2,nthroot(IRSIR_2,3));
colormap('gray');
axis square; axis xy;
xlabel('x_2 (m)'); ylabel('y_2 (m)'); title(['IRSIR_2 / I_1 = ',num2str(IRSIRvol)]);
figure(12)
plot(x_2,IRSIR_2(Ny_2/2+1,:));
xlabel('x_2 (m)'); title('IRSIR_2');



% figure(13)
% plot(x_2,unwrap(angle(URSIR_2(Ny_2/2+1,:))));
% title('wavefront');
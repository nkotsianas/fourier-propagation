% Square beam propagation example (FTFP, FIRP)

% "_1" indicates source plane variables;
% "_2" indicates observation plane variables;

wx_1 = 51e-3;
wy_1 = wx_1;
lambda = 500e-9;
z = 20000;
k = 2*pi/lambda;
lz = lambda*z;

% Here, we set Lx_1 = Lx_2 = .5;
% This results in the condition lz*Mx_2 = (Lx_2)^2;
% Given our "choices" of lambda and z, this determines Mx_2 to be 250;
% Just note that not all of these quantities are independent;
% (Similarly for the y-values);


Lx_2 = 1;
Ly_2 = 1;
Mx_2 = Lx_2^2/lz;
Ny_2 = Ly_2^2/lz;
% Mx_2 = Lx_2^2/lz + mod(Lx_2^2/lz,2);
% Ny_2 = Ly_2^2/lz + mod(Ly_2^2/lz,2);
% Mx_2 = 256;
% Ny_2 = 256;
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

Lx_1 = lz/dx_2;
Ly_1 = lz/dy_2;
Mx_1 = Mx_2;
Ny_1 = Ny_2;
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

% disp(['Mx_2 = ', num2str(Mx_2),' // ',...
%     'Lx_2 = ', num2str(Lx_2), ' // ',...
%     'Lx_1 = ', num2str(Lx_1)]);


% display some sampling parameters
disp(' wx_1 // z // Mx_1 // Lx_1');
disp([' ',num2str(wx_1),' // ',num2str(z),' // ',num2str(Mx_1),' // ',num2str(Lx_1)]);
disp('  TF<1<IR // Fraun<<Fres<1');
disp(['x: ', num2str(Lx_1^2/(Mx_1*lz)),' // ',...
    num2str(wx_1^2/lz)]);
disp(['y: ', num2str(Ly_1^2/(Ny_1*lz)),' // ',...
    num2str(wy_1^2/lz)]);
% return;
% development of irradiance plot
U_1 = Rect(X_1/(2*wx_1)).*Rect(Y_1/(2*wy_1));    % source field

figure(1)
imagesc(x_1,y_1,U_1);
colormap('gray');
axis square; axis xy;
xlabel('x_1 (m)'); ylabel('y_1 (m)'); title('U_1');
% return;




U_2 = FraunProp(U_1,Lx_1,Ly_1,lambda,z);
UTF_2 = FTFP(U_1,Lx_2,Ly_2,lambda,z);
UIR_2 = FIRP(U_1,Lx_2,Ly_2,lambda,z);
I_2 = abs(U_2).^2;
ITF_2 = abs(UTF_2).^2;
IIR_2 = abs(UIR_2).^2;


URSIR_2 = RSIRProp(U_1,Lx_1,Ly_1,lambda,z);
URSTF_2 = RSTFProp(U_1,Lx_1,Ly_1,lambda,z);
IRSIR_2 = abs(URSIR_2).^2;
IRSTF_2 = abs(URSTF_2).^2;



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





% figure(6)   % obs field x-axis phase profile
% plot(x_2,unwrap(angle(UTF_2(Ny_2/2+1,:))));
% xlabel('x_2 (m)'); title(['Phase(UTF_2) (rad); z = ',num2str(z),' m']);
% figure(7)   % obs field x-axis phase profile
% plot(x_2,unwrap(angle(UIR_2(Ny_2/2+1,:))));
% xlabel('x_2 (m)'); title(['Phase(UIR_2) (rad); z = ',num2str(z),' m']);
% 




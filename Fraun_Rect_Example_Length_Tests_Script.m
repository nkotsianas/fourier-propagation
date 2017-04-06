% Fraun_Rect_Example_Length_Tests_Script

Lxs = .1:.1:2;

outsEq = zeros(1,length(Lxs));
outsDbl = zeros(1,length(Lxs));

for k = 1:length(Lxs)
    outsEq(k) = Fraun_Rect_Example_Length_Tests(Lxs(k),Lxs(k));
    outsDbl(k) = Fraun_Rect_Example_Length_Tests(Lxs(k),2*Lxs(k));
end

figure(99)
plot(Lxs,abs(outsEq-1)+1,'-o',Lxs,abs(outsDbl-1)+1,'-o',Lxs,ones(1,length(Lxs)));
xlabel('Lx'); ylabel('I2max / I2Amax');
legend('Ly=Lx','Ly=2Lx');
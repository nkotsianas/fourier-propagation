n = 1.5;
f = 10;

dx = .01;
w = 10;
x = -w:dx:w;

s = f/(n+1)*(1 - sqrt(1 + (n+1)/(n-1)*x.^2/f^2));
sp = -s;

R = (n-1)*f;
k = -2.25;
z = x.^2./(R*(1 + sqrt(1 - (1+k)*x.^2/R^2)));
perf = f/(n-1)*(x.^2/f^2)./(1 + sqrt(1 + (n+1)/(n-1)*(x.^2/f^2)));

% asymp = f/(n+1) - x/sqrt(n^2-1);

figure(1)
plot(x,sp,x,z-.1,x,perf-.2);
return;

figure(2)
plot(x,abs(sp-z));
title(['|max| = ',num2str(max(abs(sp-z)))]);
% axis equal;
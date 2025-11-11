gamma=-20:0.01:-1;
xi=0.8;
f=((1+gamma)./gamma).^(xi+1);
g=1-(xi+1)./(xi-gamma);

figure(1)
plot(gamma,[f;g])
legend('small','big')


y=0:0.01:1;
f_bigPenis=1-y.^(xi+1);
f_small=(1-y)*(xi+1)./(xi*(1-y)+1);
figure(2)
plot(y,f_bigPenis-f_small) % must be positive

g_big=y.^(-xi);
g_small=xi*(1-y)+1;
figure(3)
plot(y,g_big-g_small)


A=2;
C=0:0.01:5;
alpha=1;
beta=2;
kappa=0.05;
l=kappa*exp((C/A).^(alpha)-(C/A).^(beta));
figure(4)
plot(C,l)

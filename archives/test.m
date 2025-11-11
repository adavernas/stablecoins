clearvars

%% PARAMETERS
r   =  0.06;
xmu =  0.05;
sig =  0.1;

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 1;
lambda = 0.02;
muk = 0.055;

varphi0 = 1;

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

%% END
writepar(mfilename)
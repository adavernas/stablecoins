

r   = 0.06;
xmu = 0.05;
sig = 0.1;

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;

xi = 10;
lambda = 0.10;

nd = 50;
dvec = linspace(0,10,nd);
gk1vec = NaN(1,nd);
gk2vec = NaN(1,nd);
gk3vec = NaN(1,nd);
fvec = NaN(1,nd);

for id=1:nd
    delt_ = dvec(id);
    akk = - sig^2/2;
    bkk =    xmu - delt_     + sig^2/2*(xi-1);
    ckk = - (xmu - delt_)*xi + sig^2/2 *xi + r + lambda - delt_ ;
    dkk = - (r   - delt_)*xi;
    
    t = (-1+sqrt(-3))/2;
    kk = [0 1 2];
    
    D0 = bkk.^2 - 3.*akk.*ckk;
    D1 = 2.*bkk.^3 - 9.*akk.*bkk.*ckk + 27.*akk.^2.*dkk;
    C = ((D1 + sqrt(D1.^2 - 4.*D0.^3))/2).^(1/3);
    
    gk = -1./(3*akk).*(bkk + t.^kk.*C + D0./(t.^kk.*C));
    
%     gam = sort(roots([akk bkk ckk dkk]),'descend');

    gk1vec(id) = gk(1);
    gk2vec(id) = gk(2);
    gk3vec(id) = gk(3);
    
    fvec(id) = sig^2/2*(xi-1) + xmu - delt_ - sig^2/2*(gk(1)+gk(2)+gk(3));
end

figure(1); clf(1); hold on;
plot(dvec,gk1vec,'m')
plot(dvec,gk2vec,'b')
plot(dvec,gk3vec,'r')

figure(2); clf(2); hold on;
plot(dvec,fvec,'k')


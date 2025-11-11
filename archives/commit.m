clearvars
set(0,'DefaultFigureWindowStyle','docked')
dbstop if error

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

% ell0 = 0.20;
% ell1 = 0.10;

% ell  = @(a) max(0,ell0 - ell1./a);
% ellp  = @(a) ell1./a.^2;
% amax  = 2*ell1/ell0;

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 6;
lambda = 0.1;

muk = 0.055;

delt0 = (r + sig^2/2 - exp(-1)*r)*1.15;
varphi0 = 0.99;

a0 = - sig^2/2;
b0 = xmu - sig^2/2;
c0 = r;

gam0 = sort(roots([a0 b0 c0]),'descend');

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

gmin = -1;
% ppmax = 400;

tolalb = 1e-6;
tolESa = 1e-5;
tolaub = 1e-4;
odeTol = 1e-8;
display = 'off';

%% END
writepar(mfilename)
par = parfun;

options = optimset('tolF',1e-20, 'tolX',1e-20);

[x1,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit fullcollateral'), amax,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit fullcollateral'), options);
[~,sol1] = Fsyst(x1,par,...
    'case','jump commit fullcollateral');
disp(sum(abs(sol1.F)))

figure(1); clf(1);
sgtitle('commit $\varphi=1$','Interpreter','LateX') 
plotfig(sol1,par)

options = optimset('tolF',1e-4,'tolX',1e-4);

LB = 1;
UB = 3;

x0    = 1.2;
alb0  = 0.02;
aub0  = 1;
EpSa0 = @(a) ones(size(a));
 
delete(gcp('nocreate'))
numCores = feature('numcores');
p = parpool(numCores);

nv = 4*numCores;
vvec = linspace(0.01,0.99,nv);
avec = linspace(1,2,30);

fvec = NaN(1,nv);
for ii=1:nv
    optClass = CommitCollateralClass();
%     optClass.albInit  = alb0;
%     optClass.aubInit  = aub0;
%     optClass.EpSaInit = EpSa0;

    fvec(ii) = optClass.optFunVarphi(x,par,vvec(ii)),x0,...
        [],[],[],[],LB,UB,[],options);
    fvec(ii) = -xF;
    avec(ii) =  x2;
    
    parsave(['./files/tmpsolution_commit_collateral_varphi_',...
        num2str(round(vvec(ii),2)*100,'%.0f'),'_' ,...
        num2str(round(avec(ii),2)*100,'%.0f')],-xF,x2)
end

kk = 0;
while 1
    kk = kk+1;
    if ~isfile(['./files/endsolution_commit_collateral_varphi_',num2str(kk),'.mat'])
        save(['./files/endsolution_commit_collateral_varphi_',num2str(kk)])
        break
    end
end

x0 = [amax;amax/2;0];

x = fsolve(@(x) Fsyst(x,par,'case','jump commit'), x0, options);
[~,sol0] = Fsyst(x,par,'case','jump commit');
disp(sum(abs(sol0.F)))
 
figure(1); clf(1);
sgtitle('commit $\varphi=0$','Interpreter','LateX') 
plotfig(sol0,par)

[x1,~] = fmincon(@(x) -FsystCon(x,par,...
    'case','jump commit fullcollateral'), amax,...
    [],[],[],[],[],[],@(x) mycon(x,par,...
    'case','jump commit fullcollateral'), options);
[~,sol1] = Fsyst(x1,par,...
    'case','jump commit fullcollateral');
disp(sum(abs(sol1.F)))

figure(2); clf(2);
sgtitle('commit $\varphi=1$','Interpreter','LateX') 
plotfig(sol1,par)

figure(3); clf(3); hold on;
plot(vvec,fvec)
plot(vvec,ones(size(vvec))*sol0.fstar)
plot(vvec,ones(size(vvec))*sol1.fstar)

title('$f^\star$','Interpreter','LateX')

figure(4); clf(4); hold on;
plot(vvec,avec)
title('$a^\star$','Interpreter','LateX')




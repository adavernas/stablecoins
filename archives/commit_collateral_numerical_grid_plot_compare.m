clearvars
set(0,'DefaultFigureWindowStyle','docked')
dbstop if error

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

savegraph = 'off';

ell1 = 1;
ell0 = ell1*r;

ell  = @(a) exp(-ell1./a).*ell0/ell1;
ellp = @(a) 1./a.^2.*exp(-ell1./a).*ell0;
amax  = ell1;

xi = 6;
lambda = 0.1;

muk = 0.055;

delt0 = (r + sig^2/2 - exp(-1)*r)*1.15;
varphi0 = 1;

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

%% END

writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20);

for hh=1:6
    figure(hh); clf(hh);
end

nl = 20;
lvec = linspace(0.01,0.5,nl);
lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
lvec = [lvec(1:34) 0.25:0.05:1];
nl = length(lvec);

nv = 20;
vvec = linspace(0.01,0.99,nv);
vvec = vvec(1:2:nv);
nv = length(vvec);

na = 10;
% avec = [linspace(1,2,na) 2.20:0.20:4];
avec = [linspace(1,2,na)];
na = length(avec);

avec_ = linspace(0,2,200);

il = 19;
iv = 3;
ia = 3;

f1 = load(['./files/savefile_function_',...
    num2str(round(avec(ia),2)*1e2,'%.0f'),...
    '_',num2str(round(vvec(iv),2)*1e2,'%.0f'),...
    '_',num2str(round(lvec(il),3)*1e3,'%.0f')]);

figure(1); hold on;
plot(avec_,f1.efun(avec_),'b')
plot(f1.aub,f1.efun(f1.aub),'*k')
plot(f1.astar,f1.efun(f1.astar),'*r')
plot(f1.astar_,f1.efun(f1.astar_),'or')
title('$e(a)$','Interpreter','LateX')

figure(2); hold on;
plot(avec_,f1.pfun(avec_),'b')
plot(f1.aub,f1.pfun(f1.aub),'*k')
plot(f1.astar,f1.pfun(f1.astar),'*r')
plot(f1.astar_,f1.pfun(f1.astar_),'or')
title('$p(a)$','Interpreter','LateX')

figure(3); hold on;
plot(avec_,f1.EpSa(avec_),'b')
plot(f1.aub,f1.EpSa(f1.aub),'*k')
plot(f1.astar,f1.EpSa(f1.astar),'*r')
plot(f1.astar_,f1.EpSa(f1.astar_),'or')
ylim([0.0 1])
title('$E[p(Sa)]$','Interpreter','LateX')

aub = fsolve(@(a) f1.eb(a,a,f1.astar),f1.aub,f1.options);
if or(aub<0,f1.eb(aub,aub,f1.astar)>1e-2)
    aub = fsolve(@(a) f1.eb(a,a,f1.astar),1,f1.options);
end

figure(4); hold on;
plot(avec_,abs(f1.eb(avec_,avec_,f1.astar)),'b')
plot(f1.aub,f1.eb(f1.aub,f1.aub,f1.astar),'*k')
title('$eb(a)$','Interpreter','LateX')

% disp(['E[p(Sastar)]: ',num2str(f1.EpSa(f1.astar_)),...
%     ' lambda: ',num2str(f1.par.lambda),...
%     ' varphi: ',num2str(f1.par.varphi0),...
%     ' fstar: ',num2str(f1.fstar),...
%     ' astar: ',num2str(f1.astar),...
%     ' error: ',num2str(f1.aout(end)-f1.aub)])
%
% il = 18;
% iv = 2;
% ia = 2;
%
% f2 = load(['./files/savefile_function_',...
%         num2str(round(avec(ia),2)*1e2,'%.0f'),...
%     '_',num2str(round(vvec(iv),2)*1e2,'%.0f'),...
%     '_',num2str(round(lvec(il),3)*1e3,'%.0f')]);
%
% figure(4); hold on;
% plot(avec_,f2.efun(avec_),'b')
% plot(f2.aub,f2.efun(f2.aub),'*k')
% plot(f2.astar,f2.efun(f2.astar),'*r')
% plot(f2.astar_,f2.efun(f2.astar_),'or')
% title('$e(a)$','Interpreter','LateX')
%
% figure(5); hold on;
% plot(avec_,f2.pfun(avec_),'b')
% plot(f2.aub,f2.pfun(f2.aub),'*k')
% plot(f2.astar,f2.pfun(f2.astar),'*r')
% plot(f2.astar_,f2.pfun(f2.astar_),'or')
% title('$p(a)$','Interpreter','LateX')
%
% figure(6); hold on;
% plot(avec_,f2.EpSa(avec_),'b')
% plot(f2.aub,f2.EpSa(f2.aub),'*k')
% plot(f2.astar,f2.EpSa(f2.astar),'*r')
% plot(f2.astar_,f2.EpSa(f2.astar_),'or')
% ylim([0.0 1])
% title('$E[p(Sa)]$','Interpreter','LateX')
%
% disp(['E[p(Sastar)]: ',num2str(f2.EpSa(f2.astar_)),...
%     ' lambda: ',num2str(f2.par.lambda),...
%     ' varphi: ',num2str(f2.par.varphi0),...
%     ' fstar: ',num2str(f2.fstar),...
%     ' astar: ',num2str(f2.astar),...
%     ' error: ',num2str(f2.aout(end)-f2.aub)])

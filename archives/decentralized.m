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
varphi0 = 0.5;

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

tolalb = 1e-6;
tolESa = 1e-5;
tolaub = 1e-4;
odeTol = 1e-7;
display = 'off';

%% END
writepar(mfilename)
par = parfun;

s = @(a) muk*ones(size(a));
delt = @(a) r-ell(a);

aa = 14;
bb = 3.5;
fontsize = 14;
titlesize = 24;
legendsize = 12;
linewidth = 1.5;

avec = linspace(0,2,200);
figure(1); clf(1); hold on;
plot(avec,s(avec),'b')
plot(avec,delt(avec),'r')

Cvec = linspace(0,2,200);
figure(2); clf(2); hold on;
h1 = plot(Cvec,s(1./Cvec),'b','linewidth',linewidth);
h2 = plot(Cvec,delt(1./Cvec),'r','linewidth',linewidth);

Cstar = fminsearch(@(C) -C*(s(1/C)-delt(1/C)),1);

x = [0           Cstar       Cstar          0];
y = [s(1./Cstar) s(1./Cstar) delt(1./Cstar) delt(1./Cstar)];
patch(x,y,[.75 0.75 0.75])

plot([Cstar Cstar],[0 delt(1./Cstar)],'--k')

xx = [0 Cstar 1 2];
set(gca, 'XTick', xx, 'xticklabel', {0, '$C^\star$', 1, 2},'TickLabelInterpreter','Latex')

legend([h1,h2],{'$s(A/C)$','$\delta(A/C)$'},...
    'Location','SouthEast','Interpreter','Latex','FontSize',legendsize)
legend boxoff

set(gca,'FontSize',fontsize)

set(gcf,'Units','inches');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition', [0 0 aa bb]);
set(gcf,'PaperSize', [aa bb]);

if ispc
    graphpath = 'E:\SHoF Dropbox\Adrien d''Avernas\Apps\Overleaf\Stable Coin\figures\'; 
elseif ismac
    graphpath = '/Users/adrien/SHoF Dropbox/Adrien d''Avernas/Apps/Overleaf/Stable Coin/figures/'; 
end

if strcmp(savegraph,'on')
    saveas(gcf,[graphpath,'decentralized'],'pdf')
end


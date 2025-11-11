clearvars
set(0,'DefaultFigureWindowStyle','docked')
addpath('./files/')
dbstop if error

%% PARAMETERS
r   = 0.06;
xmu = 0.05;
sig = 0.1;

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

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

gmin = -1;

odeTol = 1e-8;
tolalb = 1e-6;
tolESa = 1e-3;
tolaub = 1e-2;
display_plot = 'off';
display_iter = 'off';

%% END
writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20);

delete(gcp('nocreate'))

% vector for lambda
nl = 20;
lvec = linspace(0.01,0.5,nl);
lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
lvec = [lvec(1:34) 0.25:0.05:1];
nl = length(lvec);

% vector for varphi
nv_tmp = 20;
vvec_tmp = linspace(0.01,0.99,nv_tmp);
vvec_low  = vvec_tmp(1:2:nv_tmp);
vvec_high = vvec_tmp(8:nv_tmp);

% vector for alpha
na = 10;
avec_tmp = linspace(1,2,na);
avec_low = (avec_tmp(2:end)+avec_tmp(1:end-1))/2;
avec_low = avec_low(avec_low<=1.5);
avec = sort([avec_tmp avec_low]);
na = length(avec);

for il=1:nl

    % initial values
    par.lambda = lvec(il);
    x1 = fsolve(@(x) Fsyst(x, par, 'case', 'jump commit'), [amax;amax/2;0], options);
    [~,sol1] = Fsyst(x1,par,'case','jump commit');
    disp(sum(abs(sol1.F)))
    
    alb00  = 0.02;
    aub00  = sol1.aub;
    EpSa00 = sol1.EpSa;
    
    if lvec(il)<=0.15
        vvec = vvec_low;
        nv = length(vvec);
    else
        vvec = vvec_high;
        nv = length(vvec);
    end
    
    parfor iv=1:nv
        optClass = CommitCollateralClass();
        optClass.alb00 = alb00;
        optClass.aub00 = aub00;
        optClass.EpSa00 = EpSa00;
        
        for ia=1:na
            if exist(['./files/savefile_function_',...
                    num2str(round(avec(ia),2)*1e2,'%.0f'),'_',...
                    num2str(round(vvec(iv),2)*1e2,'%.0f'),'_',...
                    num2str(round(lvec(il),3)*1e3,'%.0f'),'.mat'],'file')
                
                f  = load(['./files/savefile_function_',...
                    num2str(round(avec(ia),2)*1e2,'%.0f'),'_',...
                    num2str(round(vvec(iv),2)*1e2,'%.0f'),'_',...
                    num2str(round(lvec(il),3)*1e3,'%.0f')]);
                
                f.fstar_ = (f.estar(f.aub,f.astar_)+1-f.par.varphi0)/f.astar_;
                if or(or(~strcmp(...
                        num2str(round(lvec(il),3)*1e3,'%.0f'),...
                        num2str(round(f.par.lambda,3)*1e3,'%.0f')),...
                        abs(f.aout(end)-f.aub0)/f.aub0>tolaub),...
                        abs(f.aub-f.aub0)/abs(f.aub0)>tolaub)
                    
                    optClass.optInit(lvec,il,vvec,iv,avec,ia);
                    optClass.optFunLambda(avec(ia),par,lvec,il,vvec,iv);
                end
            end
            if ~exist(['./files/savefile_function_',...
                    num2str(round(avec(ia),2)*1e2,'%.0f'),'_',...
                    num2str(round(vvec(iv),2)*1e2,'%.0f'),'_',...
                    num2str(round(lvec(il),3)*1e3,'%.0f'),'.mat'],'file')
                
                optClass.optInit(lvec,il,vvec,iv,avec,ia);
                optClass.optFunLambda(avec(ia),par,lvec,il,vvec,iv);
            end
        end
    end
end

eval('save ''./files/savefile_grid''')


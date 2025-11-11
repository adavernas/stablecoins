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
varphi0 = 1;

if r - xmu + lambda/(xi+1)<0
    disp('phi is negative!')
end

if r + xmu - sig^2 - lambda/(xi-1)<0
    disp('kappa is negative!')
end

loadsolution = 'on';
savegraph = 'off';

%% END

writepar(mfilename)
par = parfun;
options = optimset('tolF',1e-20, 'tolX',1e-20);

display = 'off';

for hh=1:10
    figure(hh); clf(hh);
end

if strcmp(loadsolution,'off')
    
    nl = 20;
    lvec = linspace(0.01,0.5,nl);
    lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
    lvec = sort([lvec (lvec(2:end)+lvec(1:end-1))/2]);
    lvec = [lvec(1:34) 0.25:0.05:1];
    nl = length(lvec);
    
    nv_tmp = 20;
    vvec_tmp = linspace(0.01,0.99,nv_tmp);
    vvec_low  = vvec_tmp(1:2:nv_tmp);
    vvec_high = vvec_tmp(8:nv_tmp);
    
    na = 10;
    % avec = [linspace(1,2,na) 2.20:0.20:4];
    avec_tmp = linspace(1,2,na);
    avec_low = (avec_tmp(2:end)+avec_tmp(1:end-1))/2;
    avec_low = avec_low(avec_low<=1.5);
    avec = sort([avec_tmp avec_low]);
    na = length(avec);
    
    x1 = [amax;amax/2;0];
    
    il0 = find(abs(lvec-0.10)==min(abs(lvec-0.10)));
    
    vstarvec_lambda = NaN(1,nl);
    fstarvec_lambda = NaN(1,nl);
    astarvec_lambda = NaN(1,nl);
    for il=1:nl
        
        if lvec(il)<=0.15
            vvec = vvec_low;
            nv = length(vvec);
        else
            vvec = vvec_high;
            nv = length(vvec);
        end
        
        astarvec = NaN(1,nv);
        fstarvec = NaN(1,nv);
        
        fvec = NaN(nv,na);
        svec = NaN(nv,na);
        uvec = NaN(nv,na);
        for iv=1:nv
            for ia=1:na
                if exist(['./files/savefile_function_',...
                        num2str(round(avec(ia),2)*1e2,'%.0f'),...
                        '_',num2str(round(vvec(iv),2)*1e2,'%.0f'),...
                        '_',num2str(round(lvec(il),3)*1e3,'%.0f'),'.mat'])
                    
                    f = load(['./files/savefile_function_',...
                        num2str(round(avec(ia),2)*1e2,'%.0f'),...
                        '_',num2str(round(vvec(iv),2)*1e2,'%.0f'),...
                        '_',num2str(round(lvec(il),3)*1e3,'%.0f')]);
                    
                    aub = fsolve(@(a) f.eb(a,a,max(a,f.astar)),f.aub,f.options);
                    if or(aub<0,f.eb(aub,aub,max(aub,f.astar))>1e-2)
                        aub = fsolve(@(a) f.eb(a,a,max(a,f.astar)),1,options);
                        if or(aub<0,f.eb(aub,aub,max(aub,f.astar))>1e-2)
                            aub = fsolve(@(a) f.eb(a,a,max(a,f.astar)),0.5,options);
                            if or(aub<0,f.eb(aub,aub,max(aub,f.astar))>1e-2)
                                aub = fsolve(@(a) f.eb(a,a,max(a,f.astar)),2,options);
                            end
                        end
                    end
                    
                    if strcmp(display,'on')
                        avec_ = linspace(0,2,100);
                        figure(1); hold on;
                        plot(avec_,f.efun(avec_),'b')
                        plot(f.aub,f.efun(f.aub),'*k')
                        plot(f.astar,f.efun(f.astar),'*r')
                        plot(f.astar_,f.efun(f.astar_),'or')
                        title('$e(a)$','Interpreter','LateX')
                        
                        figure(2); hold on;
                        plot(avec_,f.pfun(avec_),'b')
                        plot(f.aub,f.pfun(f.aub),'*k')
                        plot(f.astar,f.pfun(f.astar),'*r')
                        plot(f.astar_,f.pfun(f.astar_),'or')
                        title('$p(a)$','Interpreter','LateX')
                        
                        figure(3); hold on;
                        plot(avec_,f.EpSa(avec_),'b')
                        plot(f.aub,f.EpSa(f.aub),'*k')
                        plot(f.astar,f.EpSa(f.astar),'*r')
                        plot(f.astar_,f.EpSa(f.astar_),'or')
                        ylim([0.0 1])
                        title('$E[p(Sa)]$','Interpreter','LateX')
                        
                        figure(4); hold on;
                        plot(avec_,abs(f.eb(avec_,avec_,max(avec_,f.astar))),'b')
                        plot(f.aub,f.eb(f.aub,f.aub,max(f.aub,f.astar)),'*k')
                        title('$eb(a)$','Interpreter','LateX')
                        
                        f.fstar_ = (f.estar(f.aub,f.astar_)+1-f.par.varphi0)/f.astar_;
                        disp(['E[p(Sastar)]: ',num2str(f.EpSa(f.astar_)),...
                            ' lambda: ',num2str(f.par.lambda),...
                            ' varphi: ',num2str(f.par.varphi0),...
                            ' fstar: ',num2str(f.fstar),...
                            ' astar: ',num2str(f.astar),...
                            ' aub: ',num2str(f.aub),...
                            ' tolaub: ',num2str(f.par.tolaub),...
                            ' err aub0: ',num2str(f.aout(end)-f.aub0),...
                            ' err aub: ',num2str(abs(f.fstar_-f.fstar0)/f.fstar0)])
                    end
                    
                    if and(abs(f.aub-aub)<1e-2,...
                            and(f.fstar>0,...
                            and(abs(f.aout(end)-f.aub0)<1e-2,...
                            strcmp(num2str(round(lvec(il),3)*1e3,'%.0f'),...
                            num2str(round(f.par.lambda,3)*1e3,'%.0f')))))
                        
                        fvec(iv,ia) = f.fstar;
                        svec(iv,ia) = f.astar;
                        uvec(iv,ia) = f.aub;
                    end
                end
            end
            
            ind0 = and(fvec(iv,:)>0,and(~isnan(fvec(iv,:)),svec(iv,:)>=uvec(iv,:)));

            ifrst = min(find(ind0==1,1,'first'),find(ind0==0,1,'last'));
            if isempty(ifrst)
                ifrst = 1;
            end
            
            ilast = find(ind0==1,1,'last');
            
            if sum(ind0==1)>4
                p = polyfit(avec(ind0),fvec(iv,ind0),3);
                
                ffun = @(a) polyval(p,a);
                %             ffun = @(a) interp1(avec(ind0),fvec(iv,ind0),a,'spline');
                
                astarvec(iv) = fminbnd(@(a) -ffun(a),1,avec(ilast));
                fstarvec(iv) = ffun(astarvec(iv));
                
                if strcmp(display,'on')
                    avec_ = linspace(min(avec(ind0)),max(avec(ind0)),200);
                    figure(5); hold on;
                    plot(avec,fvec(iv,:),'ob')
                    plot(avec(ind0),fvec(iv,ind0),'ok')
                    plot(avec_,ffun(avec_))
                    plot(astarvec(iv),fstarvec(iv),'*r')
                    %                 ylim([max(max(fvec))*0.85 max(max(fvec))*1.05])
                    
                    title('$f^\star(a^\star;\varphi,\lambda)$','Interpreter','LateX')
                end
            end
            
        end
        
        iid = 1;
        if il==24
            iid=4;
        elseif il==27
            iid = 3;
        end
        
        ind1 = and(fstarvec>0,~isnan(fstarvec));
        ind1(1:iid) = 0;
        
        if sum(ind1==1)>4
            
            p = polyfit(vvec(ind1),fstarvec(ind1),3);
            ffun = @(a) polyval(p,a);
            
            p = polyfit(vvec(ind1),astarvec(ind1),3);
            afun = @(a) polyval(p,a);
            
            %         ffun = @(a) interp1(vvec(ind1),fstarvec(ind1),a,'spline');
            %         afun = @(a) interp1(vvec(ind1),astarvec(ind1),a,'spline');
            
            vstarvec_lambda(il) = fminbnd(@(a) -ffun(a),vvec(iid),vvec(end));
            fstarvec_lambda(il) = ffun(vstarvec_lambda(il));
            astarvec_lambda(il) = afun(vstarvec_lambda(il));
            
            vvec_ = linspace(0,1,100);
            
            figure(6); hold on;
            plot(vvec_,ffun(vvec_))
            plot(vvec,fstarvec,'ok')
            %     plot(1,sol0.fstar,'*k')
            %     plot(0,sol1.fstar,'*k')
            plot(vstarvec_lambda(il),fstarvec_lambda(il),'*r')
            title('$f^\star(\varphi;\lambda)$','Interpreter','LateX')
            ylim([0 2])
            
            figure(7); hold on;
            plot(vvec,astarvec)
            plot(vvec,astarvec,'ok')
            %     plot(1,sol0.astar,'*k')
            %     plot(0,sol1.astar,'*k')
            title('$a^\star(\varphi;\lambda)$','Interpreter','LateX')
            ylim([0 4])
            
            if il==il0
                ffun0 = ffun;
                fstarvec0 = fstarvec;
                vvec0 = vvec;
            end
        end
    end
    
    save('commit_collateral_solution')
else

    load('commit_collateral_solution')
end

if sum(~isnan(fstarvec_lambda))>4
    ind = find(astarvec_lambda>0);
    
    figure(7); clf(7); hold on;
    plot(lvec(ind),astarvec_lambda(ind),'ok')
    title('$a^\star(\lambda)$','Interpreter','LateX')
    
    lvec_ = linspace(0,lvec(end),200);
    
    ind = find(vstarvec_lambda>0.02);
    ind = [ind(1)-1 ind(1:end-1)];
    
    % p = polyfit([lvec(ind) lvec(end)],[vstarvec_lambda(ind) 1],2);
    % vstarfun = @(a) max(0,polyval(p, a));
    
    vstarfun = @(a) (a>=lvec(ind(1))).*interp1([lvec(ind) 1],[vstarvec_lambda(ind) 1],a,'pchip');
    
    lvec_tmp = linspace(lvec(ind(1)),1,10);
    
    vstarfun = @(a) (a>=lvec(ind(1))).*interp1(lvec_tmp,vstarfun(lvec_tmp),a,'pchip');
    
    figure(8); clf(8); hold on;
    plot(lvec_,vstarfun(lvec_))
    plot(lvec,vstarvec_lambda,'ok')
    title('$\varphi^\star(\lambda)$','Interpreter','LateX')
    ylim([0 1])
    
    sol.vstarfun = vstarfun;
    
    x0 = [amax;amax/2;0];
    
    par.lambda = 0;
    [x0,~] = fmincon(@(x) -FsystCon(x,par,...
        'case','jump commit'), x0,...
        [],[],[],[],[],[],@(x) mycon(x,par,...
        'case','jump commit'), options);
    [~,sol0] = Fsyst(x0,par,...
        'case','jump commit');
    disp(sum(abs(sol0.F)))
    
    lvec = unique(sort([lvec 0.2:0.025:0.30]));
    lvec = lvec(lvec<=0.30);
    nl = length(lvec);
    ustarvec_lambda0 = NaN(1,nl);
    astarvec_lambda0 = NaN(1,nl);
    fstarvec_lambda0 = zeros(1,nl);
    fstarvec_lambda1 = NaN(1,nl);
    
    for il=1:nl
        par.lambda = lvec(il);
        
        [x0,~] = fmincon(@(x) -FsystCon(x,par,...
            'case','jump commit'), x0,...
            [],[],[],[],[],[],@(x) mycon(x,par,...
            'case','jump commit'), options);
        [~,sol0] = Fsyst(x0,par,...
            'case','jump commit');
        disp(sum(abs(sol0.F)))
        
        if sum(abs(sol0.F))<1e-5
            ustarvec_lambda0(il) = sol0.aub;
            astarvec_lambda0(il) = sol0.astar;
            fstarvec_lambda0(il) = sol0.fstar;
        end
        
%         x1 = fsolve(@(x) Fsyst(x,par,'case','jump commit fullcollateral'), amax, options);
%         [~,sol1] = Fsyst(x1,par,'case','jump commit fullcollateral');
%         fstarvec_lambda1(il) = sol1.fstar;
    end
    sol0.fstarvec_lambda0 = fstarvec_lambda0;
    sol1.fstarvec_lambda1 = fstarvec_lambda1;
    
    ind = find(~isnan(fstarvec_lambda));
    p = polyfit([lvec(ind)],log([fstarvec_lambda(ind)]),3);
    fstarfun = @(a) exp(polyval(p, a));
    
    ind0 = find(fstarvec_lambda0>0);
    p = polyfit(lvec(ind0),fstarvec_lambda0(ind0),4);
    fstarfun0 = @(a) polyval(p, a);
    
    sol0.lvec0 = fsolve(@(a) fstarfun0(a),0.2);
    
    figure(9); clf(9); hold on;
    plot(lvec_,fstarfun(lvec_))
    plot(lvec_,fstarfun0(lvec_),'--b')
    plot(lvec,fstarvec_lambda0,'b')
    plot(lvec,fstarvec_lambda1,'k')
    plot(lvec,fstarvec_lambda,'ok')
    title('$f^\star(\lambda)$','Interpreter','LateX')
    ylim([0 2])
    
    figure(10); clf(10); hold on;
    plot(lvec(1:il),astarvec_lambda0(1:il)./ustarvec_lambda0(1:il),'b')

    title('$a^\star(\lambda)$','Interpreter','LateX')
    
    ind0 = find(and(~isnan(fstarvec_lambda),~isnan(fstarvec_lambda0)));
    
    p = polyfit(lvec(ind0),log(max(0.01,fstarvec_lambda(ind0)-fstarvec_lambda0(ind0))),3);
    sol0.dfstarfun0 = @(a) exp(polyval(p, a));
    
    p = polyfit([lvec(ind) lvec(end)],log([fstarvec_lambda(ind)-fstarvec_lambda1(ind) 0.001]),3);
    sol1.dfstarfun1 = @(a) exp(polyval(p, a));
    
    p = polyfit(lvec(ind), fstarvec_lambda(ind),3);
    sol.fstarfun = @(a) polyval(p, a);
    sol.fstarvec_lambda = fstarvec_lambda;
    
    %     figure(10); clf(10); hold on;
    %     plot(lvec,dfstarfun(lvec))
    %     plot(lvec,fstarvec_lambda-fstarvec_lambda1,'ok')
    %     title('$f^\star(\lambda)$','Interpreter','LateX')
    %     ylim([0 0.5])
    
    sol.lvec = lvec;
    
    if strcmp(savegraph,'on')
        figure(1); clf(1);
        plotgraphs(par,'graphcase','commit collateral numerical',...
            'sol',sol,'sol00',sol00,'sol0',sol0,'sol1',sol1,...
            'name','commit_collateral_numerical','savegraph','on')
    end
    
    par.lambda = 0.10;
    x2 = [1;0.5;1];
    nv = 30;
    vvec_ = linspace(0,0.98,nv);
    fstarvec_varphi0 = NaN(1,nv);
    for iv=1:nv-1
        par.varphi0 = vvec_(iv);
        
        [x2,~] = fmincon(@(x) -FsystCon(x,par,...
            'case','jump collateral delta infinite'), x2,...
            [],[],[],[],[],[],@(x) mycon(x,par,...
            'case','jump collateral delta infinite'), options);
        [~,sol2] = Fsyst(x2,par,...
            'case','jump collateral delta infinite');
        
        if sum(abs(sol2.F))<1e-5
            fstarvec_varphi0(iv) = sol2.fstar;
        end
    end
    
    figure(11); clf(11); hold on;
    plot(vvec_,ffun0(vvec_))
    plot(vvec0,fstarvec0,'ok')
    plot(vvec_,fstarvec_varphi0)
    
    sol.vvec_ = vvec_;
    sol.ffun0 = ffun0;
    sol.fstarvec_varphi0 = fstarvec_varphi0;
    
    ind = ~isnan(fstarvec_varphi0);
    p = polyfit(vvec_(ind),fstarvec_varphi0(ind),3);
    ffun = @(a) polyval(p,a);
    sol.varphi = fminbnd(@(a) -ffun(a),0,1);
    sol.fstar = ffun(sol.varphi);
    
    sol.varphi0 = fminbnd(@(a) -ffun0(a),0,1);
    sol.fstar0 = ffun0(sol.varphi0);
    
    figure(2); clf(2);
    plotgraphs(par,'graphcase','no commit collateral numerical',...
        'sol',sol,...
        'name','nocommit_collateral_numerical','savegraph',savegraph)
end

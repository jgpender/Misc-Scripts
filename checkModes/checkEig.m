% clear;

% The raw data is in ../netcdfOutput, as usual, but it's convenient to
% squash the data into single files, like this:
%   ncrcat ../netcdfOutput/TS_his_* TS_his.nc
%   ncrcat ../netcdfOutput/TS_his2_* TS_his2.nc


grid = roms_get_grid('../TS_0.25.nc','TS_his.nc',0,1);
[ny,nx] = size(grid.h);

u = nc_varget('./TS_his.nc','u_eastward');
v = nc_varget('./TS_his.nc','v_northward');
rho = nc_varget('./TS_his.nc','rho');
N2_w_roms = nc_varget('./TS_his.nc','w');

% ubar = nc_varget('./TS_his.nc','ubar_eastward');
% vbar = nc_varget('./TS_his.nc','vbar_northward');
ubar = nc_varget('./TS_his.nc','ubar');
vbar = nc_varget('./TS_his.nc','vbar');

uEig =  nc_varget('./TS_his.nc','EigC_u');
vEig =  nc_varget('./TS_his.nc','EigC_v');
rhoEig =  nc_varget('./TS_his.nc','EigC_rho');

% uEig =  nc_varget('./TS_his.nc','EigC_u');
% vEig =  nc_varget('./TS_his.nc','EigC_v');
% rhoEig =  nc_varget('./TS_his.nc','EigC_rho');

% load '../../../InputFiles/BC_IC_HYCOM_2014/PsiCalcs_001/psi.mat';
% load '../../../InputFiles/BC_IC_HYCOM_2014/PsiCalcs_001/N20.mat';

%% Verify that uEig(1) = ubar_eastward

% % fig(1);clf;pcolor(grid.lon_rho,grid.lat_rho,sq(uEig(end,1,:,:)));shading flat;colorbar;title('uEig(1)')
% % % fig(2);clf;pcolor(sq(ubar(end,:,:)));shading flat
% % fig(2);clf;pcolor(grid.lon_rho,grid.lat_rho,sq(uEig(end,1,:,:))-sq(ubar(end,:,:)));shading flat;colorbar;caxis(10^-3*[-1 1]);title('uEig(1)-ubar')

% Huh, I thought it'd be better than that. Perhaps ROMS is calculating ubar
% and vbar in a somewhat different way or at a slightly different time. Fortunately
% I have the full 3D output for u_eastward and v_eastward so I can do the
% calculation here and compare to the first velocity eigencoefficient.

myUbar = 0*sq(ubar(end,:,:));
myVbar = myUbar;
for ii=1:nx; for jj=1:ny
        dz = diff(grid.z_w(:,jj,ii));
        myUbar(jj,ii)=dot(dz,u(end,:,jj,ii))/grid.h(jj,ii);
        myVbar(jj,ii)=dot(dz,v(end,:,jj,ii))/grid.h(jj,ii);
    end;end;

% Now it works.

% % fig(3);clf;pcolor(grid.lon_rho,grid.lat_rho,sq(uEig(end,1,:,:))-sq(myUbar(:,:)));shading flat;colorbar;caxis(10^-7*[-1 1]);title('Calculate ubar discrepancy "in house"')
% % fig(4);clf;pcolor(grid.lon_rho,grid.lat_rho,sq(vEig(end,1,:,:))-sq(myVbar(:,:)));shading flat;colorbar;caxis(10^-7*[-1 1]);title('Calculate vbar discrepancy "in house"')


%% Here is a question:
%   If I were to recalculate the modes on every snapshot and then find the
%   coefficients, how different would they be?

T = nc_varget('./TS_his.nc','temp');
S = nc_varget('./TS_his.nc','salt');
[nt,nz,ny,nx]=size(S);

ii=70; jj=50;
nVel = 3;
nRho = 2;

myUeig = zeros(nt,nVel); myUeig0 = myUeig;
myVeig = zeros(nt,nVel); myVeig0 = myVeig;
myReig = zeros(nt,nRho); myReig0 = myReig;
num = myReig;
denom = myReig;
wTime = zeros(nRho,nt,nz);
pTime = zeros(nVel,nt,nz);
pXTime = zeros(nVel,nt,nz);
rXTime = zeros(nVel,nt,nz);

for tt=1:nt
    if ( grid.mask_rho(jj,ii) == 1 )
        
        z_r = -flipud(grid.z_r(:,jj,ii));
        z_w = -flipud(grid.z_w(:,jj,ii));
        dz  = diff(z_w);
        
        P_r = sw_pres(abs(z_r),-43);
        P_w = sw_pres(abs(z_w),-43);
%         T_r = flipud(T(tt,:,jj,ii)');
%         S_r = flipud(S(tt,:,jj,ii)');
        rho_r = flipud(rho(tt,:,jj,ii)');
        
        N2_w = flipud(N2_w_roms(tt,:,jj,ii)');
        
%         S_w = interp1(z_r, S_r, z_w, 'linear', 'extrap');
%         T_w = interp1(z_r, T_r, z_w, 'linear', 'extrap');
        
%         [N2_r,~,p_ave_r] = sw_bfrq(S_r,  sw_temp(S_r,T_r,P_r,0) ,P_r); N2_r(N2_r<1e-8)=1e-8;
%         [N2_w,~,p_ave_w] = sw_bfrq(S_w,  sw_temp(S_w,T_w,P_w,0) ,P_w); N2_w(N2_w<1e-8)=1e-8;
        
%         N2_r = interp1(p_ave_r,N2_r,P_r,'pchip');N2_r(N2_r<1e-8)=1e-8;
%         N2_w = interp1(p_ave_w,N2_w,P_w,'pchip');N2_w(N2_w<1e-8)=1e-8;

        N2_r = interp1(z_w, N2_w, z_r, 'linear');
        
        [w, p, ce, Pout]=ROMS_dynmodes_jgp(N2_w,P_w,P_r);
        
%         p = [ones(nz,1)/1030. p];
        pX = 0*p;
        dz= -diff( flipud(grid.z_w(:,jj,ii)));
        for nn = 1:nVel;    % nVel modes are relevant
            p(:,nn) = p(:,nn)/sign(p(1,nn));
            Pnorm = p(1,nn);
            p(:,nn) = p(:,nn)/Pnorm;
            pX(:,nn) = p(:,nn) .* dz /sum( p(:,nn).^2 .* dz );          
            
            pTime(nn,tt,:) = p(:,nn);
            pXTime(nn,tt,:) = pX(:,nn);
            
            if tt==1
                pX0=pX;
            end;
            
            myUeig(tt,nn) = dot(pX(:,nn),flipud(u(tt,:,jj,ii)'));
            myVeig(tt,nn) = dot(pX(:,nn),flipud(v(tt,:,jj,ii)'));
            
            myUeig0(tt,nn) = dot(pX0(:,nn),flipud(u(tt,:,jj,ii)'));
            myVeig0(tt,nn) = dot(pX0(:,nn),flipud(v(tt,:,jj,ii)'));
            
        end;
        
        r=0*w;
        rX=r;
        for nn = 1:nRho;    % nRow modes are relevant
            w(:,nn) = w(:,nn)/sign(w(1,nn));
            wTime(nn,tt,:) = w(:,nn);
       
            r(:,nn)  = w(:,nn).*N2_r;     
%             rX(:,nn) = (w(:,nn) .* dz) /sum( w(:,nn).^2 .* N2_r .* dz);       
            rX(:,nn) = (rho_r .* w(:,nn) .* dz) /sum( w(2:end,nn).^2 .* diff(rho_r)  );
            rXTime(nn,tt,:) = rX(:,nn);
            
            
            num(tt,nn) = dot(rho_r .* w(:,nn) .* dz,rho_r);
            denom(tt,nn) = sum( w(2:end,nn).^2 .* diff(rho_r)  );
            
            
            
% scale as per rmode
%             Cnorm = max(abs(r(:,nn)));
%             Cnorm = 10^-5;
%             r(:,nn)  = r(:,nn)/Cnorm;
%             rX(:,nn) = rX(:,nn)*Cnorm;

% scale as per r transform vector            
%             Cnorm = max(abs(rX(:,nn)))*sign(rX(1,nn));
            Cnorm = 1;
            r(:,nn)  = r(:,nn)*Cnorm;
            rX(:,nn) = rX(:,nn)/Cnorm;
            
            if tt==1
                rX0=rX;
            end;
            

            myReig(tt,nn) = dot(rX(:,nn),rho_r);
%             myReig(tt,nn) = dot(rX(:,nn),flipud(rho(tt,:,jj,ii)'));
            
%             myReig0(tt,nn) = dot(rX0(:,nn),flipud(rho(tt,:,jj,ii)'));
            myReig0(tt,nn) = dot(rX0(:,nn),rho_r);
            
        end % kk
        
%         fig(90);clf;plot(rXform(:,1,jj,ii));hold on;plot(rX(:,1),'r')
%         fig(91);clf;plot(rXform(:,2,jj,ii));hold on;plot(rX(:,2),'r')
%         fig(92);clf;plot(rX(:,1))
%         fig(93);clf;plot(rX(:,2))
        
        aaa=5;
        
        

        % Now find the eigencoefficients
        
    else
        'masked out'
    end
end; done('time')


%% How much difference does it make when you recalcultate the modes each snapshot?

% fig(1);clf;
% subplot(3,3,1);plot(myUeig(:,1));title('Recalc modes')
% subplot(3,3,2);plot(uEig(:,1,jj,ii));title('t=0 only')
% subplot(3,3,3);plot(uEig(:,1,jj,ii)-myUeig(:,1));title('difference')
% subplot(3,3,4);plot(myUeig(:,2))
% subplot(3,3,5);plot(uEig(:,2,jj,ii))
% subplot(3,3,6);plot(uEig(:,2,jj,ii)-myUeig(:,2))
% subplot(3sdf,3,7);plot(myUeig(:,3))
% subplot(3,3,8);plot(uEig(:,3,jj,ii))
% subplot(3,3,9);plot(uEig(:,3,jj,ii)-myUeig(:,3))


fig(1);clf;
subplot(3,1,1);plot(myUeig(:,1));hold on;plot(myUeig0(:,1),'r');ylabel('mode 0');title('Eig-u');
subplot(3,1,2);plot(myUeig(:,2));hold on;plot(myUeig0(:,2),'r');ylabel('mode 1')
subplot(3,1,3);plot(myUeig(:,3));hold on;plot(myUeig0(:,3),'r');ylabel('mode 2');xlabel('day');

fig(2);clf;
subplot(3,1,1);plot(myVeig(:,1));hold on;plot(myVeig0(:,1),'r');ylabel('mode 0');title('Eig-v');
subplot(3,1,2);plot(myVeig(:,2));hold on;plot(myVeig0(:,2),'r');ylabel('mode 1')
subplot(3,1,3);plot(myVeig(:,3));hold on;plot(myVeig0(:,3),'r');ylabel('mode 2');xlabel('day')

fig(3);clf;
subplot(2,1,1);plot(myReig(:,1));hold on;plot(myReig0(:,1),'r');ylabel('mode 1');title('Eig-rho');
subplot(2,1,2);plot(myReig(:,2));hold on;plot(myReig0(:,2),'r');ylabel('mode 2');xlabel('day')


fig(4);
clf;
subplot(2,2,1);plot(myReig(:,1));hold on;plot(myReig0(:,1),'r');ylabel('mode 1');title('Eig-rho');
subplot(2,2,2);plot(myReig0(:,1),'r');ylabel('mode 1');title('Eig-rho');
subplot(2,2,3);plot(myReig(:,2));hold on;plot(myReig0(:,2),'r');ylabel('mode 2');xlabel('day')
subplot(2,2,4);plot(myReig0(:,2),'r');ylabel('mode 1');title('Eig-rho');

fig(5);clf
subplot(2,1,1);plot(z_r,rho_r)
subplot(2,1,2);plot(z_r(2:end),diff(rho_r))

fig(6);clf;plot(num(:,1));title('numerator')
fig(7);clf;plot(denom(:,1));title('denominator')
fig(8);clf;plot(num(:,1) ./ denom(:,1));title('ratio')

fig(16);clf;plot(num(:,2));title('numerator')
fig(17);clf;plot(denom(:,2));title('denominator')
fig(18);clf;plot(num(:,2) ./ denom(:,2));title('ratio')

% fig(2);clf;
% subplot(3,2,1);plot(pX(:,1))
% subplot(3,2,2);plot(pX0(:,1))
% subplot(3,2,3);plot(pX(:,2))
% subplot(3,2,4);plot(pX0(:,2))
% subplot(3,2,5);plot(pX(:,3))
% subplot(3,2,6);plot(pX0(:,3))



%% Try to come up with a way to improve the velocity transform vector and the 
%   density transform vector as t progresses.
%   NOTE: I am going to fixate on the transform vectors, not the modes!

% pX transform

CpX = zeros(nVel,nt);
for tt=1:nt; for nn=1:nVel
    sq( pXTime(nn,tt,:) ./ pXTime(nn,1,:) ); CpX(nn,tt) = mean(ans);
end;end
fig(2);clf;plot(CpX(:,:)')

CrX = zeros(nRho,nt);
for tt=1:nt; for nn=1:nRho
    sq( rXTime(nn,tt,:) ./ rXTime(nn,1,:) ); CrX(nn,tt) = mean(ans);
end;end
fig(3);clf;plot(CrX(:,:)')


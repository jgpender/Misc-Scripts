clear;

grid = roms_get_grid('../TS_0.25.nc','TS_his.nc',0,1);
Z_r  = -flipdim(grid.z_r,1);
Z_w  = -flipdim(grid.z_w,1);
[nz,ny,nx]=size(Z_r);

dZ   = 0*Z_r;
for myI=1:nx; for jj=1:ny
	dZ(:,jj,myI) = diff(Z_w(:,jj,myI));
end;end; %myI,jj

Rho_r     = flipdim(nc_varget('./TS_his.nc','rho'),2);
N2_w_ROMS = flipdim(nc_varget('./TS_his.nc','w')  ,2);


%%

myT=32;myI=80;myJ=60;

myZ_r = Z_r(:,myJ,myI);
myZ_w = Z_w(:,myJ,myI);

myN2_w = sq(N2_w_ROMS(myT,:,myJ,myI));

nSmooth=2;

fig(1);clf;plot(myZ_w,myN2_w);title('N2_w');xlabel('z')
fig(2);clf;plot(myZ_w(end-10:end),myN2_w(end-10:end));title('N2_w');xlabel('z')
%     hold on;plot(myZ_w(end-10:end),smooth(myN2_w(end-10:end),nSmooth),'r')
fig(3);clf;plot(myZ_w(1:15),myN2_w(1:15));title('N2_w');xlabel('z')
%     hold on;plot(myZ_w(1:15),smooth(myN2_w(1:15),nSmooth),'r')

% myN2_w(myN2_w<1e-8)=1e-8
myN2_r = interp1(myZ_w, myN2_w, myZ_r, 'linear', 'extrap');


%%

myRho_r = sq(Rho_r(myT,:,myJ,myI));
myRho_w = interp1(myZ_r, myRho_r, myZ_w, 'linear', 'extrap');

fig(4);clf;plot(myZ_w,myRho_w,'b');%hold on;plot(myZ_r,myRho_r,'r')
title('Rho anom, as per ROMS');xlabel('z')

D = 9.8 / 1000 * diff(myRho_w) ./ myN2_r;

fig(5);clf;plot(myZ_r,D);title('D');xlabel('z')

%% Try creating a rho0 by smoothing rho

myRho_w_smoo=smooth(myRho_w,2);
myRho_r_smoo = .5*(myRho_w_smoo(1:end-1)+myRho_w_smoo(2:end));

fig(10);clf;plot(myZ_w,myRho_w,'b');hold on;plot(myZ_w,myRho_w_smoo,'r')
xlabel('z');title('ROMS rho output (blue) vs a smoothed version (red')

fig(11);clf;plot(myZ_w,myRho_w-myRho_w_smoo,'b');xlabel('z')
title('ROMS rho output MINUS a smoothed version')

newD = 9.8 * diff(myRho_w-myRho_w_smoo) ./ myRho_r_smoo ./ myN2_r
fig(12);clf;plot(myZ_r,newD);title('Recalculation of the displacement');xlabel('z')
fig(13);clf;plot(myZ_r,newD);title('Recalculation of the displacement');xlabel('z');ylim([-100 1500])





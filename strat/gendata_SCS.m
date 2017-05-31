%% This is a matlab script that generates input for the MITgcm model

% philosophy: push all generic operations like writing of outfput files
% into helper subroutines. This code should contain only experiment
% specific modifications

clear;
addpath LOCAL_MITGCM_LIBRARY % make sure we have the consistent version of the MITGCM library, superceding the ~/matlab/MITGCM version
%----------------------------------------------------------------------
% some basic parameters that can allow for easy configuration
%----------------------------------------------------------------------
model_version = ''; % a flag you can set if you want

ifcartesian  = 0; % ifcartesian = 0 is spherical (lon,lat) grid 
BT=1;nconsts = [1 2 3 4]; ifplottpxo=0;  rotateit=0;MODEL.toposmoo = 2;MODEL.topopresmoo = 3;use_sam_kelly_structure=0;OB.smoo=9;OB.flux_mag=10;

if ifcartesian 
	MODEL.dX0 = 8000;MODEL.dY0 = 8000; fac = 8000/MODEL.dX0 % model resolutions
	MODEL.rotate_angle=-32; MODEL.lon0=143+1.5;MODEL.lat0=-44.4;   % angle that cartesian grid is rotated,  geographic origin of cartesian grid
else
	MODEL.dX0 = 0.125;MODEL.dY0 = 0.125;MODEL.lon0=105;MODEL.lat0=0; fac = .125/MODEL.dX0;  % resolution and origin of model grid
end	 

MODEL.Nx=round(fac*240);MODEL.Ny=round(fac*280);                    % number of gridpoints

% JGP game delZ1 to give max(MODEL.Z) ~ 9000
vstretch = 1;	MODEL.Nz=200;MODEL.delZ0=10;MODEL.delZ1=85;
% vstretch = 1;	MODEL.Nz=200;MODEL.delZ0=10;MODEL.delZ1=110;
% vstretch = 0;	MODEL.Res_Z=50;

% adjust delXfine & delYfine to fine tune number of gridpoints in stretch
% grid
%-------------------------------------------------------------------------------
MODEL.hstretch = 0    ;STRETCH.delXfine  =   8e3;STRETCH.delYfine  =   8e3;STRETCH.delXcoarse =   8e3 ;STRETCH.xsigwidth  = 2;
STRETCH.LfineX = 200e3;STRETCH.Lcoarse1X = 75e3;STRETCH.Lcoarse1Y  = 225e3 ;%MODEL.offshoresmoo = 0;

               
% set offshoresmoo to 1 to smooth topography at edges and offshore. Only
% for cartesian grid at this time
%------------------------------------------------------
MODEL.offshoresmoo=0; 
MODEL.offshoresmoo=0; MODEL.smoofac=round(2);MODEL.eastsmoo=600e3;MODEL.northsmoo=600e3;MODEL.southsmoo=100e3;

% JGP Move calculation of zgrid to this section so that MODEL.H_min can be
% calculated
% NOTE:  When using the linear stretched grid check that MODEL.H_max is consistent 
% with the maximum value in
%   MODEL.fullZ  
% The vertical spacing is set by choosing 
%   MODEL.Nz	MODEL.delZ0	MODEL.delZ
% above.  Use
%       delZ1 = (2*MODEL.H_max - MODEL.Nz*MODEL.delZ0)/MODEL.Nz


% z-grid
%-----------------------------------------

if ~vstretch
	MODEL.Nz   =ceil(MODEL.H_max/MODEL.Res_Z);
	MODEL.delZ =ones(MODEL.Nz,1)*MODEL.Res_Z;
	MODEL.Z    =cumsum(MODEL.delZ)-MODEL.delZ/2;
else
	MODEL.delZ =linspace(MODEL.delZ0,MODEL.delZ1,MODEL.Nz)';
    MODEL.Z    =cumsum(MODEL.delZ)-MODEL.delZ/2;%[min(z) round(max(z))]
end


%JGP make a vector of all the cell edge elevations and all the cell center
%   elevations
MODEL.fullZ = sort( [MODEL.Z' cumsum(MODEL.delZ') ] )';

MODEL.H_max = 10000;
% MODEL.H_min =   50;
MODEL.H_min = sum(MODEL.delZ(1:2));


% extract the regional topography that will be regridded for the model
%------------------------------------------------
MODEL.toposrcfile = '/import/c/w/jpender/dataDir/TTide/DATA/icr_topo30_merge.mat';    
MITGCM_extract_topo(100,140,-10,40,2,'/import/c/w/jpender/dataDir/TOPO/topo30.grd',MODEL.toposrcfile);done('extract topo')
% MITGCM_extract_topo(-225+360,-180+360,-70,-30,2,'/import/c/w/jpender/dataDir/TOPO/topo30.grd',MODEL.toposrcfile);done('extract topo')
%------------------------------------------------
% specify the Cartesian grid
%------------------------------------------------
if ifcartesian
	MODEL=MITGCM_regrid_topo_cartesian(MODEL,STRETCH)       ;done('regrid topo');%return
else
	MODEL=MITGCM_regrid_topo_spherical(MODEL)
end	
%return
%-----------------------------------------
% %% z-grid
% %-----------------------------------------
% 
% if ~vstretch
% 	MODEL.Nz   =ceil(MODEL.H_max/MODEL.Res_Z);
% 	MODEL.delZ =ones(MODEL.Nz,1)*MODEL.Res_Z;
% 	MODEL.Z    =cumsum(MODEL.delZ)-MODEL.delZ/2;
% else
% 	MODEL.delZ =linspace(MODEL.delZ0,MODEL.delZ1,MODEL.Nz)';
%     MODEL.Z    =cumsum(MODEL.delZ)-MODEL.delZ/2;%[min(z) round(max(z))]
% end
% 
% 
% %JGP make a vector of all the cell edge elevations and all the cell center
% %   elevations
% MODEL.fullZ = sort( [MODEL.Z' cumsum(MODEL.delZ') ] )';




%-----------------------------------------
% stratification
%-----------------------------------------

% % Find the deepest point in the ROI   - no longer used
% [dum1,ind]=max(MODEL.H(:));
% [ind1 ind2] = ind2sub(size(MODEL.H),ind);
% MODEL.H(ind1, ind2)
% MODEL.lon_strat = MODEL.lon(ind2)
% MODEL.lat_strat = MODEL.lat(ind1)

% Choose a location near Palau that is around 5500 m deep
nlon=123;nlat=162;   % these are indices in MODEL.lon and MODEL.lat
MODEL.lon(nlon)
MODEL.lat(nlat)
MODEL.H(nlat,nlon)
MODEL.lat_strat=20.;
MODEL.lon_strat=120.3;

MODEL.N2_min=10^-6.5;


MODEL.alphaT=2e-4;
MODEL=MITGCM_get_EWG_stratification_linear_EOS_T_only(MODEL);done('getting stratification')



%!!!!!!!!!! NOTE!!!!  ROMS has the z coordinate negative in ana_initial.h
%and ana_tobc.h

MODEL.Z=-MODEL.Z





%% Try fitting a second-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),2);
Tfit=10.^(p(3)+p(2)*MODEL.Z+p(1)*MODEL.Z.^2);
fig(2);clf;
subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, second order fit')

%% Try fitting a third-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),3);
Tfit=10.^(p(4)+p(3)*MODEL.Z+p(2)*MODEL.Z.^2+p(1)*MODEL.Z.^3);
fig(3);clf;
subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, third order fit')


%% Try fitting a fourth-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),4);
Tfit=10.^(p(5)+p(4)*MODEL.Z+p(3)*MODEL.Z.^2+p(2)*MODEL.Z.^3+p(1)*MODEL.Z.^4);
fig(4);clf;
subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, fourth order fit')


%% Try fitting a fifth-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),5);
Tfit=10.^(p(6)+p(5)*MODEL.Z+p(4)*MODEL.Z.^2+p(3)*MODEL.Z.^3+p(2)*MODEL.Z.^4+p(1)*MODEL.Z.^5);
fig(5);clf;
% fig(24);clf; plot(x2JGP,y2JGP);hold on
% plot(x2JGP,log10(N2_fit2_JGP),'r');title('log10(N2.ewg) vs log10(N2.fit) third order')
% 
% fig(25);clf; plot(x2JGP,10.^(y2JGP));hold on;
% plot(x2JGP,N2_fit2_JGP,'r');title('N2.ewg vs N2.fit third order')

subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, fifth order fit')


%% Try fitting a sixth-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),6);
Tfit=10.^(p(7)+p(6)*MODEL.Z+p(5)*MODEL.Z.^2+p(4)*MODEL.Z.^3+p(3)*MODEL.Z.^4+p(2)*MODEL.Z.^5+p(1)*MODEL.Z.^6);
fig(6);clf;
% fig(24);clf; plot(x2JGP,y2JGP);hold on10.^(pJGP(6)+pJGP(5)*x2JGP+pJGP(4)*x2JGP.^2+pJGP(3)*x2JGP.^3+pJGP(2)*x2JGP.^4+pJGP(1)*x2JGP.^5)
% plot(x2JGP,log10(N2_fit2_JGP),'r');title('log10(N2.ewg) vs log10(N2.fit) third order')
% 
% fig(25);clf; plot(x2JGP,10.^(y2JGP));hold on;
% plot(x2JGP,N2_fit2_JGP,'r');title('N2.ewg vs N2.fit third order')

subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, sixth order fit')


%% Try fitting a seventh-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),7);
Tfit=10.^(p(8)+p(7)*MODEL.Z+p(6)*MODEL.Z.^2+p(5)*MODEL.Z.^3+p(4)*MODEL.Z.^4+p(3)*MODEL.Z.^5+p(2)*MODEL.Z.^6+p(1)*MODEL.Z.^7);
fig(7);clf; 
% fig(24);clf; plot(x2JGP,y2JGP);hold on10.^(pJGP(6)+pJGP(5)*x2JGP+pJGP(4)*x2JGP.^2+pJGP(3)*x2JGP.^3+pJGP(2)*x2JGP.^4+pJGP(1)*x2JGP.^5)
% plot(x2JGP,log10(N2_fit2_JGP),'r');title('log10(N2.ewg) vs log10(N2.fit) third order')
% 
% fig(25);clf; plot(x2JGP,10.^(y2JGP));hold on;
% plot(x2JGP,N2_fit2_JGP,'r');title('N2.ewg vs N2.fit third order')

subplot(1,2,1);plot(Tfit,'b');hold on;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref, seventh order fit')

p

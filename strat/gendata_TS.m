%% This is a matlab script that generates input for the MITgcm model

% philosophy: push all generic operations like writing of outfput files
% into helper subroutines. This code should contain only experiment
% specific modifications

clear;warning off
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
	MODEL.dX0 = 0.125;MODEL.dY0 = 0.125;MODEL.lon0=142;MODEL.lat0=-60; fac = .125/MODEL.dX0;  % resolution and origin of model grid
end	 

MODEL.Nx=round(fac*2*120);MODEL.Ny=round(fac*2*100);                    % number of gridpoints

vstretch = 1;	MODEL.Nz=100;MODEL.delZ0=10;MODEL.delZ1=110;
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

MODEL.H_max = 6000;
% MODEL.H_min =   50;
MODEL.H_min = sum(MODEL.delZ(1:2));

% extract the regional topography that will be regridded for the model
%------------------------------------------------
MODEL.toposrcfile = '/import/c/w/jpender/dataDir/TTide/DATA/icr_topo30_merge.mat';    
MITGCM_extract_topo(-225+360,-180+360,-70,-30,2,'/import/c/w/jpender/dataDir/TOPO/topo30.grd',MODEL.toposrcfile);done('extract topo')
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
MODEL.lon_strat = -208+360;MODEL.lat_strat = -44;

MODEL.N2_min=10^-6.5;

MODEL.alphaT=2e-4;
MODEL=MITGCM_get_EWG_stratification_linear_EOS_T_only(MODEL);done('getting stratification')

%% Try fitting a second-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),2);
Tfit=10.^(p(3)+p(2)*MODEL.Z+p(1)*MODEL.Z.^2);
fig(2);
subplot(1,2,1);plot(Tfit,'b');hold;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref')

%% Try fitting a third-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),3);
Tfit=10.^(p(4)+p(3)*MODEL.Z+p(2)*MODEL.Z.^2+p(1)*MODEL.Z.^3);
fig(3);
subplot(1,2,1);plot(Tfit,'b');hold;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref')


%% Try fitting a fourth-order polynomial to the log of Tref

p=polyfit(MODEL.Z,log10(MODEL.Tref),4);
Tfit=10.^(p(5)+p(4)*MODEL.Z+p(3)*MODEL.Z.^2+p(2)*MODEL.Z.^3+p(1)*MODEL.Z.^4);
fig(4);
subplot(1,2,1);plot(Tfit,'b');hold;plot(MODEL.Tref,'r');title('Tref and Tfit')
subplot(1,2,2);plot(Tfit-MODEL.Tref);title('Tfit-Tref')


clear;


%% JGP read in ROMS output

%ROMS.File = '../netcdfOutput/TS_his2_2014_0005.nc';
ROMS.File = 'TS_his2_2014_zeta.nc';
ROMS.zeta = nc_varget(ROMS.File,'zeta');
ROMS.ubar = nc_varget(ROMS.File,'ubar');
ROMS.vbar = nc_varget(ROMS.File,'vbar');
ROMS.time = nc_varget(ROMS.File,'ocean_time');
ROMS.timeInHours = (ROMS.time - ROMS.time(1))/3600;

fid=fopen('gfile');
ROMS.gfile = ['../',fgetl(fid)];

done('loading ROMS output')


%% Read in TPXO fields

warning('off','all');   % JGP the warning messages were driving me nutz
%----------------------------------------
% load tpxo ellipses
%----------------------------------------
%Note that OTIS screws up x & y convention. x(j) is the column  x index, y(i) is the row index:
base = '/import/c/w/jpender/ROMS/OTIS_DATA/';
ufile=[base,'u_tpxo7.2'];
gfile=[base,'grid_tpxo7.2'];
mfile=[base,'Model_tpxo7.2'];

lon0 = 142;lon1 = 192;
lat0 = -60;lat1 = -30;
hskip = 8;   %!!!!!!!!!!!!!! JGP -  this was hskip = 4

% JGP this is where you get the tide data direct from TPXO.  The footprint
% is the whole earth, not just the Tasman Sea
[tpxo.lon,tpxo.lat,tpxo.uamp,tpxo.upha]=tmd_get_coeff(mfile,'u','M2');
[tpxo.lon,tpxo.lat,tpxo.vamp,tpxo.vpha]=tmd_get_coeff(mfile,'v','M2');
[tpxo.lon,tpxo.lat,tpxo.zamp,tpxo.zpha]=tmd_get_coeff(mfile,'z','M2');


% Then confine the tpxo fields to the Tasman Sea.

tmplon=nc_varget(ROMS.gfile,'lon_rho');tmplon=tmplon(1,:);
tmplat=nc_varget(ROMS.gfile,'lat_rho');tmplat=tmplat(:,1);
idx = find(tmplon>=lon0&tmplon<=lon1);
jdx = find(tmplat>=lat0&tmplat<=lat1);
 
model.lat = tmplat(jdx(1:hskip:end));
model.lon = tmplon(idx(1:hskip:end));

idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));

tpxo.lon = tpxo.lon(idxs);
tpxo.lat = tpxo.lat(jdxs);
tpxo.zamp = tpxo.zamp(jdxs,idxs);
tpxo.zpha = tpxo.zpha(jdxs,idxs);
tpxo.uamp = tpxo.uamp(jdxs,idxs)/100;
tpxo.upha = tpxo.upha(jdxs,idxs);
tpxo.vamp = tpxo.vamp(jdxs,idxs)/100;
tpxo.vpha = tpxo.vpha(jdxs,idxs);

done('loading tpxo')



%%  Read in the ncl tidal forcing file


forcing.File = '../../../InputFiles/TidesM2/TS_tides_otps.nc';
forcing.Cmax = nc_varget(forcing.File,'tide_Cmax');
forcing.Cmin = nc_varget(forcing.File,'tide_Cmin');
forcing.Cang = nc_varget(forcing.File,'tide_Cangle');
forcing.Cpha = nc_varget(forcing.File,'tide_Cphase');


done('loading ncl forcing file')

%%  Consider the eastern boundary - start with zeta


% Here is a plot of zeta vs time on this boundary
[dumt dumy dumx] = size(ROMS.zeta);
fig(1);clf;pcolor(squeeze(ROMS.zeta(:,:,dumx))');shading flat;colorbar;hold on
    title('ROMS - zeta vs time on eastern boundary');xlabel('time (hours)');ylabel('y')

% I should be able to compare this to tpxo, right?
[dumy dumx] = size(tpxo.zamp);
tpxo.zvt=zeros(length(ROMS.time),dumy);
pi=3.14159;

for tt = 1:length(ROMS.timeInHours);for yy = 1:dumy;
      tpxo.zvt(tt,yy) = tpxo.zamp(yy,dumx)*cos(2*pi/12.42*ROMS.timeInHours(tt) - pi/180*tpxo.zpha(yy,dumx));    
    end;end
fig(2);clf;pcolor(tpxo.zvt');shading flat;colorbar;hold on
    title('TPXO - zeta vs time on eastern boundary');xlabel('time (hours)');ylabel('y')


%%  Now try ubar


% Here is a plot of zeta vs time on this boundary
[dumt dumy dumx] = size(ROMS.ubar);
fig(1);clf;pcolor(squeeze(ROMS.ubar(:,:,dumx))');shading flat;colorbar;hold on
    title('ROMS - ubar vs time on eastern boundary');xlabel('time (hours)');ylabel('y');caxis([-.015 .015]);

% I should be able to compare this to tpxo, right?
[dumy dumx] = size(tpxo.uamp);
tpxo.uvt=zeros(length(ROMS.time),dumy);
pi=3.14159;

for tt = 1:length(ROMS.timeInHours);for yy = 1:dumy;
      tpxo.uvt(tt,yy) = tpxo.uamp(yy,dumx)*cos(2*pi/12.42*ROMS.timeInHours(tt) - pi/180*tpxo.upha(yy,dumx));    
    end;end
fig(2);clf;pcolor(tpxo.uvt');shading flat;colorbar;hold on
    title('TPXO - ubar vs time on eastern boundary');xlabel('time (hours)');ylabel('y');caxis([-.015 .015]);


%%  Now try vbar


% Here is a plot of zeta vs time on this boundary
[dumt dumy dumx] = size(ROMS.vbar);
fig(1);clf;pcolor(squeeze(ROMS.vbar(:,:,dumx))');shading flat;colorbar;hold on
    title('ROMS - vbar vs time on eastern boundary');xlabel('time (hours)');ylabel('y')

% I should be able to compare this to tpxo, right?
[dumy dumx] = size(tpxo.vamp);
tpxo.vvt=zeros(length(ROMS.time),dumy);
pi=3.14159;

for tt = 1:length(ROMS.timeInHours);for yy = 1:dumy;
      tpxo.vvt(tt,yy) = tpxo.vamp(yy,dumx)*cos(2*pi/12.42*ROMS.timeInHours(tt) - pi/180*tpxo.vpha(yy,dumx));    
    end;end
fig(2);clf;pcolor(tpxo.vvt');shading flat;colorbar;hold on
    title('TPXO - vbar vs time on eastern boundary');xlabel('time (hours)');ylabel('y')







%% Plot out z,u,v for ROMS and for TPXO at a particular edge point in the fucked up zone

% ROMS.zeta is 402 points wide and 242 points high
% tpxo.zamp is 196 points wide and 115 points high

iROMS=402;jROMS=121; 
fig(1);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.zeta(:,jROMS,iROMS)));ylim([-.3 .3]); hold on

phaShift = 3.14159 /4;
itpxo=196;jtpxo=57;
plot(ROMS.timeInHours/24, tpxo.zamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.zpha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('zeta(t) - ROMS and TPXO(red)');ylim([-.3 .3])


iROMS=401;jROMS=121;
fig(2);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.ubar(:,jROMS,iROMS)));hold on
    plot(ROMS.timeInHours/24, tpxo.uamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.upha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('u(t) - ROMS and TPXO(red)');
    
    

iROMS=402;jROMS=121;
fig(3);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.vbar(:,jROMS,iROMS)));hold on
       plot(ROMS.timeInHours/24, tpxo.vamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.vpha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('v(t) - ROMS and TPXO(red)');


    

%% Plot out z,u,v for ROMS and for TPXO at a particular edge point in the OK zone

% ROMS.zeta is 402 points wide and 242 points high
% tpxo.zamp is 196 points wide and 115 points high

iROMS=402;jROMS=181; 
fig(4);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.zeta(:,jROMS,iROMS)));ylim([-.3 .3]); hold on

phaShift = -.1 *3.14159 ;
itpxo=196;jtpxo=86;
plot(ROMS.timeInHours/24, tpxo.zamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.zpha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('zeta(t) - ROMS and TPXO(red)');ylim([-.4 .4])


iROMS=401;jROMS=181;
fig(5);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.ubar(:,jROMS,iROMS)));hold on
    plot(ROMS.timeInHours/24, tpxo.uamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.upha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('u(t) - ROMS and TPXO(red)');
    
    

iROMS=402;jROMS=181;
fig(6);clf;plot(ROMS.timeInHours/24,squeeze(ROMS.vbar(:,jROMS,iROMS)));hold on
       plot(ROMS.timeInHours/24, tpxo.vamp(jtpxo,itpxo)*cos(2*3.14159/12.42*ROMS.timeInHours +...
        tpxo.vpha(jtpxo,itpxo)*3.14159/180 + phaShift),'r' );title('v(t) - ROMS and TPXO(red)');




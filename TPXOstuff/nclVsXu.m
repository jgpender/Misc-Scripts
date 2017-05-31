clear;


%% JGP read in ROMS output

ROMS.File = 'TS_his2_2014_zeta.nc';
ROMS.zeta = nc_varget(ROMS.File,'zeta');
ROMS.ubar = nc_varget(ROMS.File,'ubar');
ROMS.vbar = nc_varget(ROMS.File,'vbar');
ROMS.time = nc_varget(ROMS.File,'ocean_time');
ROMS.timeInHours = (ROMS.time - ROMS.time(1))/3600;

fid=fopen('gfile');
ROMS.gfile = ['../',fgetl(fid)]


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


% Confine the tpxo fields to the Tasman Sea.

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
tpxo.vpha = tpxo.vpha(jdxs,idxs);done('loading tpxo')


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


%%  Calculate tidal ellipse by Xu's routine, and read in the ncl results

[SEMA, ECC, INC,PHA] = ap2ep(tpxo.uamp,tpxo.upha,tpxo.vamp,tpxo.vpha);


forcing.File = '../../../InputFiles/TidesM2/TS_tides_otps.nc';
forcing.Cmax = nc_varget(forcing.File,'tide_Cmax');
forcing.Cmin = nc_varget(forcing.File,'tide_Cmin');
forcing.Cang = nc_varget(forcing.File,'tide_Cangle');
forcing.Cpha = nc_varget(forcing.File,'tide_Cphase');



%%  Issue 1:  Does Xu's Cmax look like the ncl Cmax?

SEMAmod=SEMA;SEMAmod(SEMAmod==0)=-1;
Cmaxmod=forcing.Cmax;Cmaxmod(Cmaxmod==NaN)=-1;

fig(1);clf;pcolor(SEMAmod);shading flat;colorbar;title('Cmax from Xu');caxis([-1 5])
fig(2);clf;pcolor(forcing.Cmax);shading flat;colorbar;title('Cmax from tidal forcing file');caxis([-1 5])

% Looks OK


%%  Issue 2:  Does Xu's Cmin look like the ncl Cmin?

SEMI=ECC .* SEMA;

fig(1);clf;pcolor(SEMI);shading flat;colorbar;title('Cmin from Xu');caxis([0 .25])
fig(2);clf;pcolor(forcing.Cmin);shading flat;colorbar;title('Cmin from tidal forcing file');caxis([0 .25])

% Looks OK


%%  Issue 3:  Does Xu's Cangle look like the ncl Cangle?


fig(1);clf;pcolor(INC);shading flat;colorbar;title('Cangle from Xu')
fig(2);clf;pcolor(forcing.Cang);shading flat;colorbar;title('Cangle from tidal forcing file')

% Clearly not, though I can almost convince myself all the disparities are
% shifts of 180 degrees.

% Try to rescale Xu's angle

INCshift=INC;
[dumy dumx] = size(INCshift);
for ii=1:dumx; for jj=1:dumy;
    if INCshift(jj,ii)>180
        INCshift(jj,ii)=INCshift(jj,ii)-180;
    end;
end;end;
fig(3);clf;pcolor(INCshift);shading flat;colorbar;title('Cangle from Xu')


% Try to rescale ncl's angle

NCLshift=forcing.Cang;
[dumy dumx] = size(NCLshift);
for ii=1:dumx; for jj=1:dumy;
    if NCLshift(jj,ii)<0
        NCLshift(jj,ii)=NCLshift(jj,ii)+180;
    end;
end;end;
fig(4);clf;pcolor(NCLshift);shading flat;colorbar;title('Cangle from tidal forcing file')

% Well, that's encouraging.  If you confine the ellipse angle for both Xu
% and the ncl to the range 0-180 they look about the same.




%%  Issue 4:  Does Xu's Cphase look like the ncl Cphase?


fig(1);clf;pcolor(PHA);shading flat;colorbar;title('Cphase from Xu')
fig(2);clf;pcolor(forcing.Cpha);shading flat;colorbar;title('Cphase from tidal forcing file')

% Clearly not, though I can almost convince myself all the disparities are
% shifts of 180 degrees.


%%  Walk through the ncl method

% Start with Re(u), Im(u), Re(v), Im(v).  These are:

pi=3.14159;

Re_u = tpxo.uamp .* cos(tpxo.upha *pi/180.);
Im_u = -tpxo.uamp .* sin(tpxo.upha *pi/180.);
Re_v = tpxo.vamp .* cos(tpxo.vpha *pi/180.);
Im_v = -tpxo.vamp .* sin(tpxo.vpha *pi/180.);

% fig(1);clf;pcolor(Re_u);shading flat;colorbar;caxis([-.2 .2]);
% fig(2);clf;pcolor(Im_u);shading flat;colorbar;caxis([-.2 .2]);
% fig(3);clf;pcolor(Re_v);shading flat;colorbar;caxis([-.2 .2]);
% fig(4);clf;pcolor(Im_v);shading flat;colorbar;caxis([-.2 .2]);

t1p = Re_u - Im_v;
t2p = Re_v + Im_u;
t1m = Re_u + Im_v;
t2m = Re_v - Im_u;

% fig(1);clf;pcolor(t1p);shading flat;colorbar;caxis([-.2 .2]);
% fig(2);clf;pcolor(t2p);shading flat;colorbar;caxis([-.2 .2]);
% fig(3);clf;pcolor(t1m);shading flat;colorbar;caxis([-.2 .2]);
% fig(4);clf;pcolor(t2m);shading flat;colorbar;caxis([-.2 .2]);

ap = .5 * sqrt ( t1p.^2 + t2p.^2);
am = .5 * sqrt ( t1m.^2 + t2m.^2);

% fig(1);clf;pcolor(ap+am);shading flat;colorbar;caxis([0 .25]);
% fig(2);clf;pcolor(ap-am);shading flat;colorbar;caxis([0 .25]);
% fig(3);clf;pcolor(forcing.Cmax);shading flat;colorbar;title('Cmax from tidal forcing file');caxis([0 .25])
% fig(4);clf;pcolor(forcing.Cmin);shading flat;colorbar;title('Cmin from tidal forcing file');caxis([0 .25])


ep = atan2 (t2p,t1p);
em = atan2 (t2m,t1m);

fig(1);clf;pcolor(ep);shading flat;colorbar;
fig(2);clf;pcolor(em);shading flat;colorbar;

% fig(1);clf;pcolor((ep+em)/2);shading flat;colorbar;
% fig(2);clf;pcolor((em-ep)/2);shading flat;colorbar;

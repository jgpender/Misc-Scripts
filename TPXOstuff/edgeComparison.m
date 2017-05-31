clear;

load('TS_142_192_-60_-30_4_M2.mat')



%% Plot tpxo u




figure(4);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.uamp(jdxs,idxs)/1000);axis xy;caxis([0,.1]);rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.uamp(jdxs,idxs)/1000,[0:.1:1],'LineColor','Black','Showtext','on')
               title('TPXO, M2 u amplitude');colorbar


%% JGP read in the tidal forcing file

tideFile = '../../../InputFiles/TidesM2/TS_tides_otps.nc';
forcingAmp = nc_varget(tideFile,'tide_Eamp');
forcingPhase = nc_varget(tideFile,'tide_Ephase');

[nx ny] = size(forcingAmp);

% figure(6);clf;
% imagesc([1:ny],[1:nx],forcingAmp);axis xy;caxis([0,1.2]);rect;hold on
%                contour([1:ny],[1:nx],forcingAmp,[0:.1:1],'LineColor','Black','Showtext','on')
%                title('tidal forcing file - M2 amplitude');colorbar
% figure(7);clf;
% imagesc([1:ny],[1:nx],forcingPhase);axis xy;caxis([0,360]);rect;hold on
%                contour([1:ny],[1:nx],forcingPhase,[0:30:345],'LineColor','Black','Showtext','on')
%                title('tidal forcing file - M2 phase');colorbar


%% Begin Harper's tidal component analysis

warning('off','all');   % JGP the warning messages were driving me nutz
%----------------------------------------
% load tpxo ellipses
%----------------------------------------
%Note that OTIS screws up x & y convention. x(j) is the column  x index, y(i) is the row index:
base = '/import/c/w/jpender/ROMS/OTIS_DATA/';
ufile=[base,'u_tpxo7.2'];
gfile=[base,'grid_tpxo7.2'];
mfile=[base,'Model_tpxo7.2'];
% JGP M2 only at this point
consts = {'m2'};%{'m2','s2','n2','k2','o1','p1','k1'}
%
% get model amplitude

model.exp1 = ['/import/c/w/jpender/roms-kate_svn/'];
model.exp2 = 'TS_0.125/Experiments/TS_0.125_2014_001_TidesM2_30days/';figno=1;

model.file_prefix = 'TS_';
model.file_suffix = '_2014';
disp(['!head -500 ','../log |grep "Input Grid"  | awk ''{print $4}'' >! gfile '])
eval(['!head -500 ','../log |grep "Input Grid"  | awk ''{print $4}'' >! gfile '])
fid=fopen('gfile');
model.gfile = ['../',fgetl(fid)]

%%

lon0 = 142;lon1 = 192;
lat0 = -60;lat1 = -30;
hskip = 4;   %!!!!!!!!!!!!!! JGP -  this was hskip = 4

% if ~exist( [model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
%     disp(['!ncrcat -v zeta -O ',model.exp1,model.exp2,'netcdfOutput/',model.file_prefix,'his2_2014_0*.nc ',model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
%     eval(['!ncrcat -v zeta -O ',model.exp1,model.exp2,'netcdfOutput/',model.file_prefix,'his2_2014_0*.nc ',model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
% end

if ~exist( ['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
    disp(['!ncrcat -v zeta -O ','../netcdfOutput/',model.file_prefix,'his2_2014_0*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
    eval(['!ncrcat -v "zeta,ubar,vbar" -O ','../netcdfOutput/',model.file_prefix,'his2_2014_0*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
end

%%


eval(['outfile =  ''','./',model.file_prefix,num2str(lon0),'_',num2str(lon1),'_',num2str(lat0),'_',num2str(lat1),'_',num2str(hskip),'_M2.mat'''])

% JGP this is where you get the tide data direct from TPXO
if ~exist(outfile)
[tpxo.lon,tpxo.lat,tpxo.uamp,tpxo.upha]=tmd_get_coeff(mfile,'u','M2');
[tpxo.lon,tpxo.lat,tpxo.vamp,tpxo.vpha]=tmd_get_coeff(mfile,'v','M2');
[tpxo.lon,tpxo.lat,tpxo.Uamp,tpxo.Upha]=tmd_get_coeff(mfile,'U','M2');
[tpxo.lon,tpxo.lat,tpxo.Vamp,tpxo.Vpha]=tmd_get_coeff(mfile,'V','M2');
[tpxo.lon,tpxo.lat,tpxo.amp,tpxo.pha]=tmd_get_coeff(mfile,'z','M2');done('loading tpxo')

tmplon=nc_varget(model.gfile,'lon_rho');tmplon=tmplon(1,:);
tmplat=nc_varget(model.gfile,'lat_rho');tmplat=tmplat(:,1);
idx = find(tmplon>=lon0&tmplon<=lon1);
jdx = find(tmplat>=lat0&tmplat<=lat1);
 
model.lat = tmplat(jdx(1:hskip:end));
model.lon = tmplon(idx(1:hskip:end));
% zeta = nc_varget([model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'zeta',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);

% JGP expand to include velocity fields
zeta = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'zeta',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);
ubar = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'ubar',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);
vbar = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'vbar',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);

%JGP get rid of NaNs
ubar(isnan(ubar))=0;
vbar(isnan(vbar))=0;


zetasmoo=nan*zeta;
for tdx = 1:length(zeta(:,1,1))
zetasmoo(tdx,:,:) = lowpassconv(sq(zeta(tdx,:,:)),hskip*2,hskip*2,1);
end;done('done smoothing zeta')
NX = length(model.lon);
NY = length(model.lat);
% times = datenum(1900,1,1,1,0,nc_varget([model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'ocean_time'));datestr(times(1:3))
times = datenum(1900,1,1,1,0,nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'ocean_time'));datestr(times(1:3))

%JGP expand to include velocity fields
model.amp=nan*ones(NY,NX);model.pha=model.amp;
model.ubaramp=nan*ones(NY,NX);model.ubarpha=model.ubaramp;
model.vbaramp=nan*ones(NY,NX);model.vbarpha=model.vbaramp;

for jdx = 1:NY
    disp([num2str(jdx),' of ',num2str(NY)])
for idx = 1:NX;
    dat = zeta(:,jdx,idx);
    if find(1-isnan(dat))
        
    [name,freq,tidecon,    zout]    =t_tide2(zeta    (:,jdx,idx),'start',times(1)+-1.5/24);
    [name,freq,tideconubar,ubarout] =t_tide2(ubar    (:,jdx,idx),'start',times(1)+-1.5/24);
    [name,freq,tideconvbar,vbarout] =t_tide2(vbar    (:,jdx,idx),'start',times(1)+-1.5/24);
    
% JGP      [name,freq,tidecon,zout] =t_tide2(zetasmoo(:,jdx,idx),'start',times(1)+-1.5/24);
     % find the correct tide
bingo=[];
for tdx = 1:length(name);
         tmpname = name(tdx,:);%disp([num2str(tdx) tmpname])
         if findstr(tmpname,'M2')
             bingo=tdx;
         end
end
     model.amp(jdx,idx)=tidecon(bingo,1);
     model.pha(jdx,idx)=tidecon(bingo,3);
     model.ubaramp(jdx,idx)=tideconubar(bingo,1);
     model.ubarpha(jdx,idx)=tideconubar(bingo,3);
     model.vbaramp(jdx,idx)=tideconvbar(bingo,1);
     model.vbarpha(jdx,idx)=tideconvbar(bingo,3);
    end
end
end
 eval(['save ',outfile,' model tpxo'])
 disp(['save ',outfile,' model tpxo'])
else
 disp(['loading ',outfile])
load(outfile)
end


%% zeta plots

figure(2);clf;
imagesc(model.lon,model.lat,model.amp);axis xy;caxis([0,1.2]);rect;hold on
               contour(model.lon,model.lat,model.amp,[0:.1:1],'LineColor','Black','Showtext','on')
               title('ROMS - M2 amplitude');colorbar
figure(3);clf;
imagesc(model.lon,model.lat,model.pha);axis xy;caxis([0,360]);rect;hold on
               contour(model.lon,model.lat,model.pha,[0:30:345],'LineColor','Black','Showtext','on')
               title('ROMS - M2 phase');colorbar

              

figure(4);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.amp(jdxs,idxs));axis xy;caxis([0,1.2]);rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.amp(jdxs,idxs),[0:.1:1],'LineColor','Black','Showtext','on')
               title('TPXO, M2 amplitude');colorbar
figure(5);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.pha(jdxs,idxs));axis xy;caxis([0,360]);rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.pha(jdxs,idxs),[0:30:360],'LineColor','Black','Showtext','on')
               title('TPXO, M2 phase');colorbar

%% u plots



figure(2);clf;
imagesc(model.lon,model.lat,model.ubaramp);axis xy;caxis([0,.2]);rect;hold on
               contour(model.lon,model.lat,model.ubaramp,[0:.1:1],'LineColor','Black','Showtext','on')
               title('ROMS M2 - ubar amplitude');colorbar
figure(3);clf;
imagesc(model.lon,model.lat,model.ubarpha);axis xy;caxis([0,360]);rect;hold on
               contour(model.lon,model.lat,model.ubarpha,[0:45:345],'LineColor','Black','Showtext','on')
               title('ROMS M2 - ubar phase');colorbar

figure(4);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.uamp(jdxs,idxs)/1000);caxis([0,.2]);axis xy;rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.uamp(jdxs,idxs),[0:.5:1],'LineColor','Black','Showtext','on')
               title('TPXO u, M2 amplitude');colorbar
figure(5);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.upha(jdxs,idxs));axis xy;caxis([0,360]);rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.upha(jdxs,idxs),[0:45:360],'LineColor','Black','Showtext','on')
               title('TPXO u, M2 phase');colorbar
               
%% U plots



figure(2);clf;
imagesc(model.lon,model.lat,model.ubaramp);axis xy;caxis([0,.2]);rect;hold on
               contour(model.lon,model.lat,model.ubaramp,[0:.1:1],'LineColor','Black','Showtext','on')
               title('ROMS M2 - ubar amplitude');colorbar
figure(3);clf;
imagesc(model.lon,model.lat,model.ubarpha);axis xy;caxis([0,360]);rect;hold on
               contour(model.lon,model.lat,model.ubarpha,[0:45:345],'LineColor','Black','Showtext','on')
               title('ROMS M2 - ubar phase');colorbar

figure(4);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.Uamp(jdxs,idxs)/1000);caxis([0,.2]);axis xy;rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.Uamp(jdxs,idxs),[0:.5:1],'LineColor','Black','Showtext','on')
               title('TPXO U, M2 amplitude');colorbar
figure(5);clf;
idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.Upha(jdxs,idxs));axis xy;caxis([0,360]);rect;hold on
               contour(tpxo.lon(idxs),tpxo.lat(jdxs),tpxo.Upha(jdxs,idxs),[0:45:360],'LineColor','Black','Showtext','on')
               title('TPXO U, M2 phase');colorbar
               
               
%%
return
   

%%

% 
% forcingFile = '/import/c/w/jpender/ROMS/realWork/EIO_32Grid/Tides_M2only_Sindhu/EIO_32_tides_otps_topo30.nc';
% forcingAmp = nc_varget(forcingFile,'tide_Eamp');
% figure(2);clf;
% hformat(12);
% 
% imagesc(forcingAmp);axis xy;caxis([0,1.8]);rect;hold on
%                contour(forcingAmp,0:.1:1,'LineColor','Black','Showtext','on');hold on
%                contour(forcingAmp,1.1:.1:1.8,'LineColor','Black')
%                title('forcing amplitude, M2 forcing');colorbar




%%
f3;clf;wysiwyg; lon0 = 92;lat0 = 5;ampmax=2;
subplot(2,2,1); % 
 idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end))
 jdxs = nearest(tpxo.lat,lat0);
 plot(tpxo.lon(idxs),tpxo.pha(jdxs(1),idxs),'r');ylim([0,360]);hold on
  jdxs = nearest(model.lat,lat0);
 plot(lons,model.pha(jdxs,:));hold on;
 title('phase along lat0')
subplot(2,2,2) % 
 idxs = nearest(tpxo.lon,lon0);
 jdxs = find(tpxo.lat>=0&tpxo.lat<=24);
 plot(tpxo.lat(jdxs),tpxo.pha(jdxs,idxs),'r');ylim([0,360]);hold on
 idxs = nearest(lons,lon0);
 jdxs = find(model.lat>=0&model.lat<=24);
 plot(model.lat(jdxs),model.pha(jdxs,idxs),'b');ylim([0,360]);hold on
 title('phase along lon0');axis tight;
subplot(2,2,3)% 
 idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end))
 jdxs = find(tpxo.lat>=lat0&tpxo.lat<=lat0)
 plot(tpxo.lon(idxs),tpxo.amp(jdxs(1),idxs),'r');hold on;;ylim([0,1])%ylim([0,360])
 jdxs = nearest(model.lat,lat0);
 plot(lons,model.amp(jdxs,:));hold on;
 title('amp along lat0')
subplot(2,2,4) ;
 idxs = nearest(tpxo.lon,lon0);
 jdxs = find(tpxo.lat>=0&tpxo.lat<=24);
 plot(tpxo.lat(jdxs),tpxo.amp(jdxs,idxs),'r');hold on;ylim([0,1])
 idxs = nearest(lons,lon0);
 jdxs = find(model.lat>=0&lats<=24);
 plot(lats(jdxs),model.amp(jdxs,idxs),'b');;ylim([0,1])
 title('amp along lon0');axis tight;ylim([0,ampmax])
%%
%%
                                               tide     freq       amp     amp_err    pha    pha_err     snr
%t_tide(zetasmoo(:,jdx,idx) )                  *M2   0.0805114    0.3793    0.017   185.09     2.73  4.9e+02
%t_tide(zetasmoo(:,jdx,idx),'start',times(1))  *M2   0.0805114    0.3793    0.017    23.69     2.44  5.1e+02
%t_tide(zetasmoo(:,jdx,idx),'start',times(1),'rayleigh',['M2  ']) M2   0.0805114    0.3793    0.004    23.68     0.52  9.3e+03
%%
%%
clear ss
tfile = '/import/c/w/hsimmons/DATA/TOPO/topo30.nc';
hskip=8;
tmplon = nc_varget(tfile,'lon');
tmplat = nc_varget(tfile,'lat');
idx=find(tmplon>=model.lon(1)&tmplon<=model.lon(end));
jdx=find(tmplat>=model.lat(1)&tmplat<=model.lat(end));
ss.lon=nc_varget(tfile,'lon',idx(1)-1,floor(length(idx)/hskip),hskip);
ss.lat=nc_varget(tfile,'lat',jdx(1)-1,floor(length(jdx)/hskip),hskip);
ss.z  =nc_varget(tfile,'z',[jdx(1)-1,idx(1)-1],floor([length(jdx),length(idx)]/hskip),[hskip hskip]);
ss.z(ss.z>0)=0;ss.z=abs(ss.z);done
%%
f1;clf;
%imagesc(model.lon,model.lat,model.h);rect;axis xy;caxis
contour(   ss.lon,   ss.lat,lowpass2d(   ss.z,3,3),[1 1],'g','linew',2);rect;axis xy;caxis;hold on
contour(model.lon,model.lat,lowpass2d(model.h,3,3),[20 30 40],'b');rect;axis xy;caxis;hold on
contour(   ss.lon,   ss.lat,lowpass2d(   ss.z,3,3),[20 30 40],'k');rect;axis xy;caxis;hold on;grid on
%%
%%
f2;clf;%wysiwyg
%imagesc(model.lon,model.lat,model.h);rect;axis xy;caxis
subplot(2,2,1);imagesc(model.lon,model.lat,lowpass2d(model.h,1,1));hold on
               contour(   ss.lon,   ss.lat,lowpass2d(   ss.z,1,1),[1 1],'g','linew',2);rect;axis xy;caxis;hold on
               hlim(95,98,13.5,17.5);caxis([0,200]);colorbar;title('model');vlines([97])
subplot(2,2,2);imagesc(   ss.lon,   ss.lat,lowpass2d(   ss.z,1,1));hold on
               contour(   ss.lon,   ss.lat,lowpass2d(   ss.z,1,1),[1 1],'g','linew',2);rect;axis xy;caxis;hold on
               hlim(95,98,13.5,17.5);caxis([0,200]);colorbar;title('smith and sandwell');;vlines([97])
subplot(2,1,2)
lon0= 97;idxm = nearest(model.lon,lon0);idxt = nearest(   ss.lon,lon0);
 plot(model.lat,model.h(:,idxm),'r');hold on
 plot(   ss.lat,   ss.z(:,idxt),'b');hold on;axis ij tight;legend('model','smith and sandwell',4)
 xlim([10.5,17.5]);ylim([0,200]);xlabel('lat');ylabel('depth');title('depth along 97E')
%% way too slow.
% lons = lon0:lon1;
% lats = lat0:lat1;
% TS=ones(length(times),length(lats),length(lons));
% for jdx =1:length(lats);for idx=1:length(lons);
%  [TS(:,jdx,idx),ConList]=tmd_tide_pred(mfile,times,lats(jdx),lons(idx),'z',1);
% end;end
%%
for us = {'u','v'}
  for contno = 1:length(consts)
  tmpu = char(us);  tmpconst=char(consts(contno));
  disp(['[tpxo_lon,tpxo_lat,amp',tmpu,',phase',tmpu,']=tmd_get_coeff(mfile,''',tmpu,''',''',tmpconst,''');'])       
  eval(['[tpxo_lon,tpxo_lat,amp',tmpu,',phase',tmpu,']=tmd_get_coeff(mfile,''',tmpu,''',''',tmpconst,''');'])       
  eval(['amp'  ,tmpu,'_',tmpconst,' = amp'  ,tmpu,'(:,:);'])
  eval(['phase',tmpu,'_',tmpconst,' = phase',tmpu,'(:,:);'])
  clear ampv ampu phaseu phasev
  end
end
%%
% get TPXO bathy
[tmp,tmp,tpxoD]=get_grid(mfile);tpxoD=tpxoD';

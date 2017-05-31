clear;
warning('off','all');   % JGP the warning messages were driving me nutz
%----------------------------------------
% load tpxo ellipses
%----------------------------------------
%Note that OTIS screws up x & y convention. x(j) is the column  x index, y(i) is the row index:
base = '/import/c1/w/jpender/ROMS/OTIS_DATA/';
ufile=[base,'u_tpxo7.2'];
gfile=[base,'grid_tpxo7.2'];
mfile=[base,'Model_tpxo7.2'];

[~,dum] = unix('ls ../*.nc');  model.gfile = dum(1:end-1)
[~,dum] =  unix('ls ../netcdf_All/*_his_* | head -1'); model.HISfile = dum(1:end-1)

grid = roms_get_grid(model.gfile,model.HISfile,0,1);

%% Here is the stuff you want to edit

[~,dum] = unix('ls ../netcdf_All | tail -1 | cut -d "_" -f1');model.file_prefix = [dum(1:end-1),'_']
model.file_suffix = '';

% PALAU_120B footprint
lon0 = 131.5;lon1 = 136.5;
lat0 = 5.5;  lat1 = 10;
hskip = 4;   %!!!!!!!!!!!!!! JGP -  this was hskip = 4

myTide = 'M2';
% Single tidal component
consts = {lower(myTide)};%{'m2','s2','n2','k2','o1','p1','k1'}


% make a custom mask_rho based on a new D_min
D_min   = 1000;
model.h = nc_varget(model.gfile,'h');
mask_rho = nc_varget(model.gfile,'mask_rho');
% fig(1);clf;pcolor(model.h);shading flat;colorbar
% fig(2);clf;pcolor(mask_rho);shading flat;colorbar
mask_rho(model.h < D_min) = 0; mask_rho(mask_rho == 0) = nan;mask_rho=mask_rho(2:end-1,2:end-1);
% fig(3);clf;pcolor(mask_rho);shading flat;colorbar


%% Archive zeta, ubar and vbar in their own file if not already done.
if ~exist( ['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
    disp(['!ncrcat -v zeta -O ','../netcdf_All/',model.file_prefix,'his2_*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
    eval(['!ncrcat -v zeta -O ','../netcdf_All/',model.file_prefix,'his2_*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'])
    done('writing zeta')
    eval(['!ncrcat -v ubar -O ','../netcdf_All/',model.file_prefix,'his2_*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_ubar.nc'])
    done('writing ubar')
    eval(['!ncrcat -v vbar -O ','../netcdf_All/',model.file_prefix,'his2_*.nc ','./',model.file_prefix,'his2',model.file_suffix,'_vbar.nc'])
    done('writing vbar')
end

eval(['outfile =  ''','./',model.file_prefix,num2str(lon0),'_',num2str(lon1),'_',num2str(lat0),'_',num2str(lat1),'_',num2str(hskip),'_',myTide,'.mat'''])

%% Do the tidal analysis (if it's not already done)

if ~exist(outfile)
    
    % TPXO stuff
    [tpxo.lon,tpxo.lat,tpxo.amp_zeta,tpxo.pha_zeta]=tmd_get_coeff(mfile,'z',myTide);
    [tpxo.lon,tpxo.lat,tpxo.amp_ubar,tpxo.pha_ubar]=tmd_get_coeff(mfile,'u',myTide);
    [tpxo.lon,tpxo.lat,tpxo.amp_vbar,tpxo.pha_vbar]=tmd_get_coeff(mfile,'v',myTide);done('loading tpxo')
    
    % ROMS stuff
    
    tmplon=nc_varget(model.gfile,'lon_rho');tmplon=tmplon(1,:);
    tmplat=nc_varget(model.gfile,'lat_rho');tmplat=tmplat(:,1);
    idx = find(tmplon>=lon0&tmplon<=lon1); model.idx = idx;
    jdx = find(tmplat>=lat0&tmplat<=lat1); model.jdx = jdx;
    
    
    % Get the area differentials for Harper's figure of merit.
    %   Try to do this right. Even if lat and lon are on a perfectly
    %   rectilinear grid, this means dx and dy will vary. If it's a
    %   telescoping grid then this is all the more true.
    %   roms_get_grid has given me x_rho, y_psi, etc etc
    
    X  = grid.x_rho(1:4:end,1:4:end);
    Y  = grid.y_rho(1:4:end,1:4:end);
    dX = X(:,2:end) - X(:,1:end-1);
    dY = Y(2:end,:) - Y(1:end-1,:);
    model.dA = dX(1:end-1,:) .* dY(:,1:end-1);
    
    model.lat = tmplat(jdx(1:hskip:end));
    model.lon = tmplon(idx(1:hskip:end));
    model.mask_rho  = mask_rho(1:hskip:end,1:hskip:end);
    
    % Interpolate the relevant part of the TPXO data onto the ROMS xy grid.
    % Also, convert from cm/s to m/s.
    
    idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end)); tpxo.idxs = [min(idxs)-1,idxs,max(idxs)+1]; %tpxo.idxs=idxs;
    jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end)); tpxo.jdxs = [min(jdxs)-1,jdxs,max(jdxs)+1]; %tpxo.jdxs=jdxs;
    
    tpxo.amp_zeta_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.amp_zeta(tpxo.jdxs,tpxo.idxs),model.lon,model.lat);
    tpxo.pha_zeta_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.pha_zeta(tpxo.jdxs,tpxo.idxs),model.lon,model.lat);
    tpxo.amp_ubar_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.amp_ubar(tpxo.jdxs,tpxo.idxs),model.lon,model.lat)/100;
    tpxo.pha_ubar_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.pha_ubar(tpxo.jdxs,tpxo.idxs),model.lon,model.lat)/100;
    tpxo.amp_vbar_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.amp_vbar(tpxo.jdxs,tpxo.idxs),model.lon,model.lat)/100;
    tpxo.pha_vbar_BB = interp2(tpxo.lon(tpxo.idxs),tpxo.lat(tpxo.jdxs),tpxo.pha_vbar(tpxo.jdxs,tpxo.idxs),model.lon,model.lat)/100;
    
    % Load in zeta, then filter with the enlarged land mask
    % Note that Harper isn't loading in the very edge of zeta so it's sized
    % smaller by 2 elements in each direction.
    zeta = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'zeta',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);
    % zetasmoo=nan*zeta;
    [nt,ny,nx] = size(zeta);
    for tt = 1:nt
        zeta(tt,:,:) = sq(zeta(tt,:,:)) .* model.mask_rho;
    end;
    done('zeta')
    
    % Load in ubar and vbar, then filter with the enlarged land mask. Ubar is on the
    % ugrid which is a leetle bit bigger than the rho grid, but if I use
    % Harper's limits in the nc_varget call I'll end up with a properly sized
    % array. Same with vbar and the v grid.
    ubar = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_ubar.nc'],'ubar',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);
    vbar = nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_vbar.nc'],'vbar',[0,jdx(1)-1,idx(1)-1],ceil([-1,length(model.lat),length(model.lon)]),[1,hskip,hskip]);
    for tt = 1:nt
        ubar(tt,:,:) = sq(ubar(tt,:,:)) .* model.mask_rho;
        vbar(tt,:,:) = sq(vbar(tt,:,:)) .* model.mask_rho;
    end;
    done('ubar and vbar')
    
    
    
    
    
    % for tdx = 1:length(zeta(:,1,1))
    % % zetasmoo(tdx,:,:) = lowpassconv(sq(zeta(tdx,:,:)),hskip*2,hskip*2,1);
    % end;done('done smoothing zeta')
    NX = length(model.lon);
    NY = length(model.lat);
    % times = datenum(1900,1,1,1,0,nc_varget([model.exp1,model.exp2,'TPXOstuff/',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'ocean_time'));datestr(times(1:3))
    times = datenum(1900,1,1,1,0,nc_varget(['./',model.file_prefix,'his2',model.file_suffix,'_zeta.nc'],'ocean_time'));datestr(times(1:3))
    
    
    
    % zeta section
    
    model.amp_zeta=nan*ones(NY,NX);model.pha_zeta=model.amp_zeta;
    for jdx = 1:NY
        disp([num2str(jdx),' of ',num2str(NY),' zeta'])
        for idx = 1:NX;
            dat = zeta(:,jdx,idx);
            if find(1-isnan(dat))
                [name,freq,tidecon,zout] =t_tide2(zeta    (:,jdx,idx),'start',times(1)+-1.5/24);
                % JGP      [name,freq,tidecon,zout] =t_tide2(zetasmoo(:,jdx,idx),'start',times(1)+-1.5/24);
                % find the correct tide
                bingo=[];
                for tdx = 1:length(name);
                    tmpname = name(tdx,:);%disp([num2str(tdx) tmpname])
                    if findstr(tmpname,myTide)
                        bingo=tdx;
                    end
                end
                model.amp_zeta(jdx,idx)=tidecon(bingo,1);
                model.pha_zeta(jdx,idx)=tidecon(bingo,3);
            end
        end
    end
    
    
    % ubar section
    
    model.amp_ubar=nan*ones(NY,NX);model.pha_ubar=model.amp_ubar;
    for jdx = 1:NY
        disp([num2str(jdx),' of ',num2str(NY),' ubar'])
        for idx = 1:NX;
            dat = ubar(:,jdx,idx);
            if find(1-isnan(dat))
                [name,freq,tidecon,zout] =t_tide2(ubar    (:,jdx,idx),'start',times(1)+-1.5/24);
                % find the correct tide
                bingo=[];
                for tdx = 1:length(name);
                    tmpname = name(tdx,:);%disp([num2str(tdx) tmpname])
                    if findstr(tmpname,myTide)
                        bingo=tdx;
                    end
                end
                model.amp_ubar(jdx,idx)=tidecon(bingo,1);
                model.pha_ubar(jdx,idx)=tidecon(bingo,3);
            end
        end
    end
    
    
    % vbar section
    
    model.amp_vbar=nan*ones(NY,NX);model.pha_vbar=model.amp_vbar;
    for jdx = 1:NY
        disp([num2str(jdx),' of ',num2str(NY),' vbar'])
        for idx = 1:NX;
            dat = vbar(:,jdx,idx);
            if find(1-isnan(dat))
                [name,freq,tidecon,zout] =t_tide2(vbar    (:,jdx,idx),'start',times(1)+-1.5/24);
                % find the correct tide
                bingo=[];
                for tdx = 1:length(name);
                    tmpname = name(tdx,:);%disp([num2str(tdx) tmpname])
                    if findstr(tmpname,myTide)
                        bingo=tdx;
                    end
                end
                model.amp_vbar(jdx,idx)=tidecon(bingo,1);
                model.pha_vbar(jdx,idx)=tidecon(bingo,3);
            end
        end
    end
    
    
    %% Harper's figure of merit
    
    % It'd be nice if I could put an image into these scripts. A picture's
    % worth a thousand words....
    % Anyhoo, it looks like the square root of
    %
    %   <  doubleIntegralOverArea(zeta_roms - zeta_tpxo)^2  >
    % The area differential
    %   dA = dx dy
    % is constant for this particular grid, but I want to write this up so that
    % it works on one of the telescoping grids. Assuming I've dont this
    % correctly, I have the area differential sized the same as zeta_roms and
    % zeta_tpxo.
    
    % size(model.amp_zeta)
    % size(tpxo.amp_zeta_BB)
    % size(model.mask_rho)
    % size(model.mask_rho)
    
    num=(model.amp_zeta - tpxo.amp_zeta_BB).^2 .*model.dA.*model.mask_rho;
    den=model.dA.*model.mask_rho;
    model.D_zeta = sqrt( nansum(num(:))/nansum(den(:)) )
    
    num=(model.amp_ubar - tpxo.amp_ubar_BB).^2 .*model.dA.*model.mask_rho;
    den=model.dA.*model.mask_rho;
    model.D_ubar = sqrt( nansum(num(:))/nansum(den(:)) )
    
    num=(model.amp_vbar - tpxo.amp_vbar_BB).^2 .*model.dA.*model.mask_rho;
    den=model.dA.*model.mask_rho;
    model.D_vbar = sqrt( nansum(num(:))/nansum(den(:)) )
    
    % sqrt(sum(ans(:))/sum(model.dA(:)))
    % fig(99);pcolor(ans);shading flat;colorbar
    
    %%
    eval(['save ',outfile,' model tpxo'])
    disp(['save ',outfile,' model tpxo'])
else
    disp(['loading ',outfile])
    load(outfile)
end



%%  JGP plots take place here

% unix('\rm -r figures');
unix('mkdir figures');


% Zeta tides

figure(1);clf;
subplot(2,2,1);  myMean = nanmean(model.amp_zeta(:));
imagesc(model.lon,model.lat,model.amp_zeta);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,model.amp_zeta,[0:.05:1],'LineColor','Black','Showtext','on')
title( ['ROMS - ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .55])


subplot(2,2,2)
imagesc(model.lon,model.lat,model.pha_zeta);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,model.pha_zeta,[0:30:345],'LineColor','Black','Showtext','on')
title(['ROMS - ',myTide,' phase']);colorbar

subplot(2,2,3);
% idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
% jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
tpxo.amp_zeta_BB; myMean=nanmean(ans(:));
imagesc(model.lon,model.lat,tpxo.amp_zeta_BB);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,tpxo.amp_zeta_BB,[0:.05:1],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .55])

subplot(2,2,4)
% idxs = find(tpxo.lon>=model.lon(1)&tpxo.lon<=model.lon(end));
% jdxs = find(tpxo.lat>=model.lat(1)&tpxo.lat<=model.lat(end));
imagesc(model.lon,model.lat,tpxo.pha_zeta_BB);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,tpxo.pha_zeta_BB,[0:30:360],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' phase']);colorbar

suptitle(['Zeta    D=',num2str(model.D_zeta)])
print(['figures/Zeta_',myTide,'tides'],'-djpeg');


%% Ubar tides

figure(2);clf;
subplot(2,2,1);  myMean = nanmean(model.amp_ubar(:));
imagesc(model.lon,model.lat, model.amp_ubar);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,model.amp_ubar,[0:.05:1],'LineColor','Black','Showtext','on')
title( ['ROMS - ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .035])

subplot(2,2,2)
imagesc(model.lon,model.lat,model.pha_ubar);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,model.pha_ubar,[0:30:345],'LineColor','Black','Showtext','on')
title(['ROMS - ',myTide,' phase']);colorbar

subplot(2,2,3); tpxo.amp_ubar_BB; myMean=nanmean(ans(:));
imagesc(model.lon,model.lat,tpxo.amp_ubar_BB);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,tpxo.amp_ubar_BB,[0:.05:1],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .035])

subplot(2,2,4)
imagesc(model.lon,model.lat,tpxo.pha_ubar_BB);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,tpxo.pha_ubar_BB,[0:30:360],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' phase']);colorbar

suptitle(['Ubar    D=',num2str(model.D_ubar)])
print(['figures/Ubar_',myTide,'tides'],'-djpeg');

%% Vbar tides

figure(3);clf;
subplot(2,2,1);  myMean =  nanmean(model.amp_vbar(:));
imagesc(model.lon,model.lat,  model.amp_vbar);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,model.amp_vbar,[0:.05:1],'LineColor','Black','Showtext','on')
title(['ROMS - ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .035])

subplot(2,2,2)
imagesc(model.lon,model.lat,model.pha_vbar);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,model.pha_vbar,[0:30:345],'LineColor','Black','Showtext','on')
title(['ROMS - ',myTide,' phase']);colorbar

subplot(2,2,3); tpxo.amp_vbar_BB; myMean=nanmean(ans(:));
imagesc(model.lon,model.lat,tpxo.amp_vbar_BB);axis xy;caxis([0,.3]);rect;hold on
contour(model.lon,model.lat,tpxo.amp_vbar_BB,[0:.05:1],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' amp (mean=',num2str( myMean ),')']);colorbar;caxis([0 .035])

subplot(2,2,4)
imagesc(model.lon,model.lat,tpxo.pha_vbar_BB);axis xy;caxis([0,360]);rect;hold on
contour(model.lon,model.lat,tpxo.pha_vbar_BB,[0:30:360],'LineColor','Black','Showtext','on')
title(['TPXO, ',myTide,' phase']);colorbar

suptitle(['Vbar    D=',num2str(model.D_vbar)])
print(['figures/Vbar_',myTide,'tides'],'-djpeg');



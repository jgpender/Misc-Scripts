clear;

%% sfos grid stuff

unix('\rm -r figures');
unix('mkdir figures');

[dum1,dum2]=unix('pwd | rev |cut -d "/" -f2 | rev | cut -d "_" -f1,2');
sfos.myGrid = dum2(1:end-1);

[dum1,dum2]=unix('ls .. |grep nc');
sfos.myGridFile = ['../',dum2(1:end-1)];

[dum1,dum2]=unix('pwd | rev |cut -d "/" -f2 | rev | cut -d "_" -f1');
sfos.myPrefix = lower(dum2(1:end-1));

HISfiles = dir('../netcdf_All/*_his_*');

sfos.grid = roms_get_grid(sfos.myGridFile,['../netcdf_All/',HISfiles(1).name],0,1);

%% Box in the voxels

% MITgcm does a better job of making the following point: the modeling volume
% consists of a bunch of building blocks all stacked up. The dimensions can
% vary from place to place. The weird thing is that there is something
% called the psi grid, and if you look at lat_psi and lon_psi these are
% exactly the lat/lon you had in mind when you created the grid in the
% first place. The other grids (lat_rho, lon_v, etc) are often offset from
% the psi grid because the points lie on one cell face or the other. In
% fact, the rho grid lies entirely on the cell faces so is larger by
% 1 point in both x and y than the number of cells.

% Most, if not all, of the variables I'm interested in are specified on the 
% rho grid. What I will do is pick a psi-grid point, average T (or
% whatever) on the 4 rho-grid faces of my cell, then assign the averaged T
% to my psi cell.

% dz is easy. Both grid.z_r and grid.z_w are on the horizontal rho grid so
% I just have to diff(grid.z_w) to get grid.dz_r. This works out because
% the length of z_w is always one greater than the length of z_r.

[nzRho,nyRho,nxRho] = size(sfos.grid.z_r);


% Here is dz on the rho grid

sfos.grid.dz_rho = zeros(nzRho,nyRho,nxRho);
for ii=1:nxRho; for jj=1:nyRho;
	sfos.dz_rho(:,jj,ii) = diff(sfos.grid.z_w(:,jj,ii));
end;end;


% Here is dz on the psi grid. Do this with a nanmean of the 4 adjacent
% values of dz on the rho grid.

sfos.dz_psi=zeros(nzRho,nyRho-1,nxRho-1);
for ii=1:nxRho-1; for jj=1:nyRho-1;
    sfos.dz_psi(:,jj,ii) = nanmean( [sfos.dz_rho(:,jj,ii),sfos.dz_rho(:,jj+1,ii),sfos.dz_rho(:,jj,ii+1),sfos.dz_rho(:,jj+1,ii+1) ],2);
end;end;



% Here is a generalized psi mask to make it easier to do the calculation.

sfos.mask_psi3D = zeros(nzRho,nyRho-1,nxRho-1);
for kk=1:nzRho;
    sfos.mask_psi3D(kk,:,:) = sfos.grid.mask_psi;
end;


% Here is a version of the rho mask that has NaNs for zeros

mask_rhoNaN = sfos.grid.mask_rho;
mask_rhoNaN(mask_rhoNaN == 0) = NaN;
sfos.mask_rho3D = zeros(nzRho,nyRho,nxRho);
for kk=1:nzRho;
    sfos.mask_rho3D(kk,:,:) = mask_rhoNaN;
end;



% Calculate delta_x for all the psi voxels. NOTE: I know that diff is
% supposed to be able to dosfos.basinArea =  this, but it's more transparent to my eye if I
% do this with a couple for loops

diff_x_rho = diff(sfos.grid.x_rho',1)';
diff_y_rho = diff(sfos.grid.y_rho,1);

dx_psi = .5 * ( diff_x_rho(1:end-1,:) + diff_x_rho(2:end,:) );
dy_psi = .5 * ( diff_y_rho(:,1:end-1) + diff_y_rho(:,2:end) );

% Turn dx and dy into wasted-space 3D matrices

sfos.dx_psi = zeros(nzRho,nyRho-1,nxRho-1);
sfos.dy_psi = zeros(nzRho,nyRho-1,nxRho-1);
for kk=1:nzRho;
    sfos.dx_psi(kk,:,:) = dx_psi;
    sfos.dy_psi(kk,:,:) = dy_psi;
end;

done('calculating dz_psi, etc')

%% Calculate basin volume and area.

% sfos

% basin area first
% 
% sfos.dx_psi .* sfos.dy_psi;sq(ans(1,:,:));
% sfos.basinAreaFull = sum(ans(:));
% inKmSq = sqrt( sfos.basinAreaFull /10^6);

% disp(['Full basin area is a square about ',num2str(inKmSq),' km on a side']);

sfos.dx_psi .* sfos.dy_psi .* sfos.mask_psi3D;sq(ans(end,:,:));
sfos.basinArea = sum(ans(:));
disp(['The water surface area is ',num2str(sfos.basinArea /10^6),' km^2']);

% basin volume next

% sfos.dx_psi .* sfos.dy_psi .* sfos.dz_psi;
% sfos.basinVolumeFull = sum(ans(:));
% disp(['Full basin volume is ',num2str(sfos.basinVolumeFull/10^9),' km^3']);

sfos.dx_psi .* sfos.dy_psi .* sfos.dz_psi .* sfos.mask_psi3D;
sfos.basinVolume = sum(ans(:));
disp(['The water volume is ',num2str(sfos.basinVolume/10^9),' km^3']);


%% Find all the averages

% Do all the snapshots or just do the first few snapshots

Nt = length(HISfiles);
% Nt = 5;

% put the state variables onto the psi grid by doing a nan over the 4
% adjacent points on the rho grid or the 2 adjacent points on the u or v
% grid.


for tt=1:Nt
    
    [dum1,dum2] = roms_get_date(['../netcdf_All/',HISfiles(tt).name]);
    sfos.time(tt)       = dum1;
    dum2
    
    % read in the fields from the netcdf file.
    
    Tdz_rho  = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'temp') )  .* sfos.mask_rho3D .* sfos.dz_rho;
    Sdz_rho  = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'salt') )  .* sfos.mask_rho3D .* sfos.dz_rho;   
  
    SSH_rho    = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'zeta') )    .* sfos.grid.mask_rho;   
    shflux_rho = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'shflux') )  .* sfos.grid.mask_rho;
    ssflux_rho = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'ssflux') )  .* sfos.grid.mask_rho;
    swrad_rho  = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'swrad') )   .* sfos.grid.mask_rho;
    
    sustr_u    = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'sustr') )   .* sfos.grid.mask_u;
    svstr_v    = sq(nc_varget(['../netcdf_All/',HISfiles(tt).name],'svstr') )   .* sfos.grid.mask_v;

    done('reading in data')
% This method does the 4-corner mean the "right way"    
%     for ii=1:nxRho-1; for jj=1:nyRho-1;
%     	SSH_psi(jj,ii)   = nanmean( [SSH_rho(  jj,ii),SSH_rho(  jj+1,ii),SSH_rho(  jj,ii+1),SSH_rho(  jj+1,ii+1) ]  );
%     	Tdz_psi(:,jj,ii) = nanmean( [Tdz_rho(:,jj,ii),Tdz_rho(:,jj+1,ii),Tdz_rho(:,jj,ii+1),Tdz_rho(:,jj+1,ii+1) ],2);
%     	Sdz_psi(:,jj,ii) = nanmean( [Sdz_rho(:,jj,ii),Sdz_rho(:,jj+1,ii),Sdz_rho(:,jj,ii+1),Sdz_rho(:,jj+1,ii+1) ],2);
%     end;end;
    
% This method is much faster but will enlarge the land mask a little. This
% only takes place in very shallow water so won't have much of an effect on
% the final answer.

    % Put the data on the psi grid.

	Tdz_psi = ( Tdz_rho(:,1:end-1,1:end-1) + Tdz_rho(:,2:end,1:end-1) + Tdz_rho(:,1:end-1,2:end) + Tdz_rho(:,2:end,2:end) ) /4;
	Sdz_psi = ( Sdz_rho(:,1:end-1,1:end-1) + Sdz_rho(:,2:end,1:end-1) + Sdz_rho(:,1:end-1,2:end) + Sdz_rho(:,2:end,2:end) ) /4;    
    
	SSH_psi    = ( SSH_rho(   1:end-1,1:end-1) + SSH_rho(   2:end,1:end-1) + SSH_rho(   1:end-1,2:end) + SSH_rho(   2:end,2:end) ) /4;
	shflux_psi = ( shflux_rho(1:end-1,1:end-1) + shflux_rho(2:end,1:end-1) + shflux_rho(1:end-1,2:end) + shflux_rho(2:end,2:end) ) /4;
	ssflux_psi = ( ssflux_rho(1:end-1,1:end-1) + ssflux_rho(2:end,1:end-1) + ssflux_rho(1:end-1,2:end) + ssflux_rho(2:end,2:end) ) /4;
	swrad_psi  = (  swrad_rho(1:end-1,1:end-1) +  swrad_rho(2:end,1:end-1) +  swrad_rho(1:end-1,2:end) +  swrad_rho(2:end,2:end) ) /4;
    
	sustr_psi  = (  sustr_u(1:end-1 ,   :     ) + sustr_u(2:end ,   :    ) ) /2;
	svstr_psi  = (  svstr_v( :      ,  1:end-1) + svstr_v( :    ,  2:end ) ) /2;    

    done('putting data on psi grid')
     
            Tdz_psi .* sfos.mask_psi3D    .*    sfos.dx_psi         .*    sfos.dy_psi        ;      
	sfos.Temp_BasinAve(tt) = nansum(ans(:)) / sfos.basinVolume;
    
            Sdz_psi .* sfos.mask_psi3D    .*    sfos.dx_psi         .*    sfos.dy_psi        ;      
	sfos.Salt_BasinAve(tt) = nansum(ans(:)) / sfos.basinVolume;
    
    
           sq(Tdz_psi(end,:,:)) ./ sq(sfos.dz_psi(end,:,:)) .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.SST_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;    
    
            sq(Sdz_psi(end,:,:)) ./ sq(sfos.dz_psi(end,:,:)) .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.SSS_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;   
    
    
    
    
    
    
            SSH_psi    .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.SSH_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;
    
        
            shflux_psi .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.shflux_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;  
    
            ssflux_psi .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.ssflux_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea; 
    
            swrad_psi  .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.swrad_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;  
    
            sustr_psi  .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.sustr_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;  
    
            svstr_psi  .* sfos.grid.mask_psi .* sq(sfos.dx_psi(1,:,:)) .* sq(sfos.dy_psi(1,:,:));      
	sfos.svstr_BasinAve(tt) = nansum(ans(:)) / sfos.basinArea;  

    
    done('finding basin averages')
    
end;


%% Save the important stuff

save('basinAve.mat','-struct','sfos')


%% Plots


fig(1);plot(sfos.time,sfos.Temp_BasinAve);title('Basin-average temperature');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/Temp_BasinAve','-djpeg');

fig(2);plot(sfos.time,sfos.Salt_BasinAve);title('Basin-average salinity');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/Salt_BasinAve','-djpeg');

fig(3);plot(sfos.time,sfos.SSH_BasinAve);title('Basin-average zeta');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/SSH_BasinAve' ,'-djpeg');

fig(4);plot(sfos.time,sfos.SST_BasinAve);title('Basin-average surface temperature');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/SST_BasinAve' ,'-djpeg');

fig(5);plot(sfos.time,sfos.SSS_BasinAve);title('Basin-average surface salinity');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/SSS_BasinAve' ,'-djpeg');



fig(6);plot(sfos.time,sfos.shflux_BasinAve);title('Basin-average shflux');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/shflux_BasinAve' ,'-djpeg');

fig(7);plot(sfos.time,sfos.ssflux_BasinAve);title('Basin-average ssflux');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/ssflux_BasinAve' ,'-djpeg');

fig(8);plot(sfos.time,sfos.swrad_BasinAve);title('Basin-average swrad');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/swrad_BasinAve' ,'-djpeg');

fig(9);plot(sfos.time,sfos.sustr_BasinAve);title('Basin-average sustr');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/sustr_BasinAve' ,'-djpeg');

fig(10);plot(sfos.time,sfos.svstr_BasinAve);title('Basin-average svstr');xlabel('date -   2016')
datetick('x','mm/dd','keepticks')
print('figures/svstr_BasinAve' ,'-djpeg');



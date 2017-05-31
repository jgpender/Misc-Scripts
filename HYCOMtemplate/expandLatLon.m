clear

%% First do the grid file

gridFile = 'HYCOM_GLBa0.08_GUAM_grid.nc';

lat = flipud(nc_varget(gridFile,'lat1D'));
lon = nc_varget(gridFile,'lon1D');
z   = nc_varget(gridFile,'z');

nx = length(lon);
ny = length(lat);
nz = length(z);

newLat = repmat(lat,1,nx);
newLon = (repmat(lon,1,ny))';

% Variables section

dum.Name = 'lon';
dum.Nctype = 'double';
dum.Dimension = {'lat','lon'};
dum.Attribute = struct('Name',{'long_name','units'},'Value',{'Longitude','degrees_east'});
nc_addvar(gridFile,dum);

dum.Name = 'lat';
dum.Nctype = 'double';
dum.Dimension = {'lat','lon'};
dum.Attribute = struct('Name',{'long_name','units'},'Value',{'Latitude','degrees_north'});
nc_addvar(gridFile,dum);

nc_varput(gridFile,'lon',newLon);
nc_varput(gridFile,'lat',newLat);

%% Now do the data files

HYCOMnames = dir('./data');

for ii=3:length(HYCOMnames)
    dataFile = ['./data/',HYCOMnames(ii).name];
    
    lat = flipud(nc_varget(dataFile,'lat1D'));
    lon = nc_varget(dataFile,'lon1D');
    z   = nc_varget(dataFile,'z');
    
    nx = length(lon);
    ny = length(lat);
    nz = length(z);
    
    newLat = repmat(lat,1,nx);
    newLon = (repmat(lon,1,ny))';
    
    % Variables section
    
    dum.Name = 'lon';
    dum.Nctype = 'double';
    dum.Dimension = {'lat','lon'};
    dum.Attribute = struct('Name',{'long_name','units'},'Value',{'Longitude','degrees_east'});    
    nc_addvar(dataFile,dum);
    
    dum.Name = 'lat';
    dum.Nctype = 'double';
    dum.Dimension = {'lat','lon'};
    dum.Attribute = struct('Name',{'long_name','units'},'Value',{'Latitude','degrees_north'});
    nc_addvar(dataFile,dum);
    
    nc_varput(dataFile,'lon',newLon);
    nc_varput(dataFile,'lat',newLat);
    
    
end;




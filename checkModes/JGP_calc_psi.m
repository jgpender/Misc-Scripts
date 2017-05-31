clear;

dirname=pwd;

[~,dum2]=unix('pwd |rev | cut -d "_" -f1 | rev');
myFiles.day = dum2(1:end-1);

[~,dum2]=unix(['ls ../ini_',myFiles.day,' |grep HY']);
myFiles.IC = ['../ini_',myFiles.day,'/',dum2(1:end-1)]
% myFiles.IC = ['TS_his_3.nc']

[~,dum2]=unix(['ls |grep TS | grep 1.nc']);
myFiles.typHIS = dum2(1:end-1);

cd ../../..
[~,dum2]=unix('pwd |rev | cut -d "/" -f1 | rev');
myFiles.expt = dum2(1:end-1);

cd(dirname);

myFiles.grid = ['../../Gridpak/',myFiles.expt,'.nc'];

myFiles

% load testVar.mat

        

%%
roms.grd = roms_get_grid(myFiles.grid,myFiles.typHIS,0,1);%done('roms_get_grid')

[~,~,roms.grd.dzr] = roms_zint(nc_varget(myFiles.typHIS,'temp'),roms.grd);

z_r = -flipdim(roms.grd.z_r,1);
z_w = -flipdim(roms.grd.z_w,1);
dzr =  flipdim(roms.grd.dzr,1);

H   = roms.grd.h;
mask= roms.grd.mask_rho;
% clearvars roms;

N20file = 'N20.mat';
myFiles.type = 'tidesAll';
psiFile = 'psi.mat';

%     if ~exist(N20file)
T = flipdim(nc_varget(myFiles.IC,'temp'),1);done('loading temp for psi')
if strmatch(myFiles.type,'tideonly')
    S = 34.7+0*T;done('simple EOS, salt = 34.7 for psi')
else
    S = flipdim(nc_varget(myFiles.IC,'salt'),1);done('loading salt for psi')
end
[nz ny nx] = size(T);

%%
if ~exist(N20file)
    
    P_r   = zeros(nz,ny,nx);
    N20_r = P_r;
    rho_r = P_r;
    
    P_w   = zeros(nz+1,ny,nx);
    N20_w = P_w;
    
    for ii = 1:nx;disp(['calculating N20 in psi, ',num2str(ii),' of ',num2str(nx)])
        for jj = 1:ny
          
          if  mask(jj,ii)==1;
            P_r(:,jj,ii) = sw_pres(abs(z_r(:,jj,ii)),-43);  % JGP calculate P_r for archive
            P_w(:,jj,ii) = sw_pres(abs(z_w(:,jj,ii)),-43);  % JGP use P_w in dynmodes
            
            tmpS_r = S(:,jj,ii);
            tmpT_r = T(:,jj,ii);
            tmpP_r = P_r(:,jj,ii);
            
            tmpz_r = z_r(:,jj,ii);
            tmpz_w = z_w(:,jj,ii);
            
            % JGP note: interp1 will give a warning message if tmpS_r has NaNs, which
            % is the case for any point under the land mask
            
            tmpS_w = interp1(tmpz_r, tmpS_r, tmpz_w, 'linear', 'extrap');
            tmpT_w = interp1(tmpz_r, tmpT_r, tmpz_w, 'linear', 'extrap');
            tmpP_w = P_w(:,jj,ii);
            
            %   fig(1);clf;plot(z(:,jj,ii),tmpS_r,'b');hold on;plot(zw(:,jj,ii),tmpS_w,'r');title('S_r vs S_w')
            %   fig(2);clf;plot(z(:,jj,ii),tmpT_r,'b');hold on;plot(zw(:,jj,ii),tmpT_w,'r');title('T_r vs T_w')
            
            if max(tmpS_r)>0;
                
                [tmpN2_r,~,p_ave_r] = sw_bfrq(tmpS_r,  sw_temp(tmpS_r,tmpT_r,tmpP_r,0) ,tmpP_r); tmpN2_r(tmpN2_r<1e-8)=1e-8;
                [tmpN2_w,~,p_ave_w] = sw_bfrq(tmpS_w,  sw_temp(tmpS_w,tmpT_w,tmpP_w,0) ,tmpP_w); tmpN2_w(tmpN2_w<1e-8)=1e-8;
                
                tmpN2_r = interp1(p_ave_r,tmpN2_r,tmpP_r,'pchip');tmpN2_r(tmpN2_r<1e-8)=1e-8;
                tmpN2_w = interp1(p_ave_w,tmpN2_w,tmpP_w,'pchip');tmpN2_w(tmpN2_w<1e-8)=1e-8;
                
                N20_r(:,jj,ii) = tmpN2_r;
                N20_w(:,jj,ii) = tmpN2_w;
            
                rho_r(:,jj,ii) = sw_dens(tmpS_r, tmpT_r, tmpP_r);
                
                %    fig(3);clf;plot(tmpP_r,tmpN2_r,'b');hold on;plot(tmpP_w,tmpN2_w,'r');title('N20_r vs N20_w')
                
            end % if mask  
                
            end % if
        end % jj
    end % ii
    eval(['save -v7.3 ',N20file,' N20_r P_r N20_w P_w rho_r'])
else
    eval(['load       ',N20file,' N20_w N20_r P_w P_r'])
end
done('N2')
%%

if ~exist(psiFile)

nm=49;

% pmodes=nan*ones(nz,nm,ny,nx);
pmodes=zeros(nz,nm,ny,nx);
rmodes=zeros(nz,nm,ny,nx);
ce    =zeros(   nm,ny,nx);
% JGP add
pXform=zeros(nz,nm,ny,nx);
rXform=zeros(nz,nm,ny,nx);

%%
%keyboard
%%

% JGP NOTE
%   T and S are on the z_rho grid (50 elements) but the wmode BC are applied
%   on the z_w grid (51 elements) at z=0 and z=H. Harper deals with this in
%   dynmodes by prepending z=0 onto the rho grid data. dynmodes returns the
%   wmodes and pmodes as 51-element vectors so the first element gets
%   stripped off.
%
%   I use interp1 to get N20 and P on the z_w grid and pass that to
%   dynmodes, which then honestly calculates the wmodes on the z_w grid (51
%   elements). diff(wmodes)/diff(z_w) gives the pmodes on the z_rho grid
%   without any fuss. interp1 puts the wmodes onto the z_rho grid.

tic
for ii = 1:nx;disp(['calculating psi ',num2str(ii),' of ',num2str(nx)])
    for jj = 1:ny
        if max(P_w(:,jj,ii))>50;
%         if find(min(N20_w(:,jj,ii)>0)&max(P_w(:,jj,ii))>50)
            
%     jj=10;ii=10;

            
            % JGP note. In MITgcm they defined the grid in terms of rectilinear
            % cells. The cell face (top and bottom) determined the w grid and the
            % cell centers determined the rho grid. Is this strictly the case with
            % ROMS? It turns out that this is NOT strictly true, thought is not
            % tooooo bad. An illustration follows:
            %     ( zw(1:5,jj,ii) + zw(2:6,jj,ii))/2   -    z(1:5,jj,ii)
            %     fig(99);clf;plot(  ( zw(1:end-1,jj,ii) + zw(2:end,jj,ii))/2   -    z(:,jj,ii)   )
            
            % JGP This returns the modes on the z_rho grid
            %   [ tmpwmodes,tmppmodes,tmpce,Pout]=ROMS_dynmodes(N20(:,jj,ii),P(:,jj,ii));
            [ tmpwmodes,tmppmodes,tmpce,Pout]=ROMS_dynmodes_jgp(N20_w(:,jj,ii),P_w(:,jj,ii),P_r(:,jj,ii));
            
            % JGP add barotropic mode to the pmodes then scale/normalize
            tmppmodes = [ones(nz,1)/1030. tmppmodes];
            pNorm = zeros(nz,1);
            for nn = 1:nm+1;
                tmppmodes(:,nn) = tmppmodes(:,nn)/sign(tmppmodes(1,nn));
                Pnorm(nn) = tmppmodes(1,nn);
                tmppmodes(:,nn) = tmppmodes(:,nn)/Pnorm(nn);
            end;     
            
            for nn = 1:nm
                tmpwmodes(:,nn) = tmpwmodes(:,nn)/sign(tmpwmodes(1,nn));
            end;

            % Find the p transform
            tmpdzr = dzr(:,jj,ii);
            for nn = 1:nm;
                pXform(:,nn,jj,ii) = tmppmodes(:,nn) .* tmpdzr /sum( tmppmodes(:,nn).^2 .* tmpdzr);              
                pmodes(:,nn,jj,ii) = tmppmodes(:,nn);
            end;
%             pXform(:,:,jj,ii)=inv(tmppmodes');
            
            % some checking
%             pinv=inv(tmppmodes);
%             fig(1);plot(pXform(:,1:5,jj,ii))
%             fig(2);plot(pinv(1:5,:)')
%             ndum=5; dum=[1:5];
%             for nn=1:ndum
%                 dum(nn) = dot(pXform(:,nn,jj,ii),u(:,jj,ii));
%             end;
%             Clin=pmodes(:,:,jj,ii)\u(:,jj,ii);
%             dum
%             Clin(1:ndum)'
%             
            aaa=5;
            
            
            % find the rho transform
            for nn = 1:nm;
                rmodes(:,nn,jj,ii)=tmpwmodes(:,nn).*N20_r(:,jj,ii);             
                rXform(:,nn,jj,ii) = (tmpwmodes(:,nn) .* tmpdzr) /sum( tmpwmodes(:,nn).^2 .* N20_r(:,jj,ii) .* tmpdzr);
                Cnorm = max(abs(rmodes(:,nn,jj,ii)));
%                 Cnorm = 1;
                rmodes(:,nn,jj,ii)=rmodes(:,nn,jj,ii)/Cnorm;
                rXform(:,nn,jj,ii) = rXform(:,nn,jj,ii)*Cnorm;
            end % kk
            
%             % some checking
%             ndum=5; dum=[1:5];
%             for nn=1:ndum
%                 dum(nn) = dot(rXform(:,nn,jj,ii),rho(:,jj,ii));
%             end;
%             Clin=rmodes(:,:,jj,ii)\rho(:,jj,ii);
%             dum
%             Clin(1:ndum)'      
            
%             dumwmodes = [ones(nz,1) tmpwmodes ];
            


% fig(1);clf;plot(-fliplr(roms.grd.z_r(:,jj,ii)'),tmppmodes(:,1:5))
% fig(2);clf;plot(-roms.grd.z_r(:,jj,ii),rmodes(:,1:5,jj,ii))   
% fig(3);clf;plot(-fliplr(roms.grd.z_r(:,jj,ii)'),tmpwmodes(:,1:5))
% fig(4);clf;plot(-fliplr(roms.grd.z_r(:,jj,ii)'),N20_r(:,jj,ii))
            aaa=5;
            

            
            ce(:,jj,ii)=tmpce(1:nm);
            
        end % N2
    end % jj%%


%!!!!!!!!!!!!! Update for each variable
end % ii
toc

%%
disp(['save -v7.3 ',psiFile,' pmodes rmodes pXform rXform ce nm'])
eval(['save -v7.3 ',psiFile,' pmodes rmodes pXform rXform ce nm'])


else
    eval(['load       ',psiFile,' pmodes rmodes pXform rXform ce nm'])
end
done('psi')
%%

nVel=3;

for ii=1:nVel
    JGP_writeNC( ['XFMvel',num2str(ii)], sq(pXform(:,ii,:,:)) )
end

nRho=2;

for ii=1:nRho
    JGP_writeNC( ['XFMrho',num2str(ii)], sq(rXform(:,ii,:,:)) )
end


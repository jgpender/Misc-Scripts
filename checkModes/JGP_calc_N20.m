function JGP_calc_N20(z_r,z_w,S,T,mask,N20file)

[nz ny nx] = size(z_r);

P_r   = nan*ones(nz,ny,nx);
N20_r = nan*ones(nz,ny,nx);
rho_r = nan*ones(nz,ny,nx);

P_w   = nan*ones(nz+1,ny,nx);
N20_w = nan*ones(nz+1,ny,nx);

for ii = 1:nx;disp(['calculating N20 in psi, ',num2str(ii),' of ',num2str(nx)])
    for jj = 1:ny
        if ( mask(jj,ii) == 1 )
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
            
            rho_r(:,jj,ii) = sw_dens(tmpS_r, tmpT_r, tmpP_r)
            
            
            %    fig(3);clf;plot(tmpP_r,tmpN2_r,'b');hold on;plot(tmpP_w,tmpN2_w,'r');title('N20_r vs N20_w')
            
        end % if
        
        end % if
        
    end % jj
end % ii
eval(['save -v7.3 ',N20file,' N20_r P_r N20_w P_w rho_r'])


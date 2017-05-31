function [wmodes_r,pmodes_r,ce,pout]=dynmodes_jgp(Nsq,p_w,p_r)
% DYNMODES calculates ocean dynamic vertical modes
%  taking a column vector of Brunt-Vaisala values (Nsq) at
%  different pressures (p) and calculating some number of 
%  dynamic modes (nmodes). 
%  Note: The input pressures need not be uniformly spaced, 
%    and the deepest pressure is assumed to be the bottom.
%
%  USAGE: [wmodes,pmodes,ce]=dynmodes(Nsq,p,nmodes);
%                               or
%                            dynmodes;  % to show demo 
%
%     Inputs: 	Nsq = column vector of Brunt-Vaisala buoyancy frequency (s^-2)
%		    	  p = column vector of pressure (decibars)
%           nmodes = number of vertical modes to calculate 
%  
%       Outputs:   wmodes = vertical velocity structure
%                  pmodes = horizontal velocity structure
%                      ce = modal speed (m/s)
%  developed by J. Klinck. July, 1999
%  send comments and corrections to klinck@ccpo.odu.edu


% JGP Modifications
% PRESUME N2 and P are on the w grid

% JGP diagnostic section as per HLS
% % % if nargin<1
% % %    help(mfilename);
% % %    nplot=3;
% % % %    test problems
% % % %      problem 1
% % % %    solution is h = ho sin(z /ce) where ce = 1 / n pi
% % % %     ce = 0.3183 / n = [ 0.3183 0.1591 0.1061]
% % % %p=0:.05:1;
% % % %z=-p;
% % % %n=length(p);
% % % %Nsq(1:n)=1;
% % % %
% % % %      problem 2
% % % %    solution is h = ho sin(No z /ce) where ce = No H / n pi
% % % %    for No=1.e-3 and H = 400, the test values are 
% % % %     ce = 0.127324 / n = [ 0.127324, 0.063662, 0.042441]
% % % %
% % % 	p=0:10:400;
% % % 	z=-p;
% % % 	n=length(p);
% % % 	Nsq(1:n)=1.e-6;
% % % 
% % % 	nmodes=3;
% % % 
% % % 	[wmodes,pmodes,ce]=dynmodes(Nsq,p,nmodes);
% % % 
% % % 	figure(1)
% % % 	plot(Nsq,z);
% % % 	title('Buoyancy Frequency Squared (s^{-2})')
% % % 
% % % 	figure(2)
% % % 	plot(ce(1:nplot),'r:o');
% % % 	title(' Modal Speed (m/s)')
% % % 
% % % 	figure(3)
% % % 	plot(wmodes(:,1:nplot),z);
% % % 	title('Vertical Velocity Structure')
% % % 
% % % 	figure(4)
% % % 	plot(pmodes(:,1:nplot),z);
% % % 	title('Horizontal Velocity Structure')
% % % 
% % %         figure(gcf)
% % %         return
% % % end
  
rho0=1028;

%    convert to column vector if necessary
[m,n] = size(p_w);
if n == 1
   p_w=p_w';
end
[m,n] = size(p_r);
if n == 1
   p_r=p_r';
end
[m,n] = size(Nsq);
if n == 1
   Nsq=Nsq';
   n=m;
end

% JGP should be using w-grid, so p(1)=0
%                 check for surface value
% % % if p_w(1) > 0
% % % %             add surface pressure with top Nsq value
% % %     z(1)=0;
% % %     z(2:n+1)=-p_w(1:n);
% % %     N2(1)=Nsq(1);
% % %     N2(2:n+1)=Nsq(1:n);
% % %     nz=n+1;
% % % else
% % %     z=-p_w;
% % %     N2=Nsq;
% % %     nz=n;
% % % end

z_w=-p_w;       % 51 elements
z_r=-p_r;       % 50 elements
% z_w=abs(p_w);       % 51 elements
% z_r=abs(p_r);       % 50 elements
% N2_w=Nsq;       % 51 elements
N2_w=Nsq;       % 51 elements
N2_w(1)=N2_w(2);
nz_w=n;         % nz_w=51
nz_r=n-1;       % nz_w=50

%  Just use z_w and z_r

% % % dz(1:nz-1)=z(1:nz-1)-z(2:nz);
% % % %        midpoint depth
% % % zm(1:nz-1)=z(1:nz-1)-.5*dz(1:nz-1)'';
% % % %        midpoint spacing
% % % dzm=zeros(1,nz);
% % % dzm(2:nz-1)=zm(1:nz-2)-zm(2:nz-1);
% % % dzm(1)=dzm(2);
% % % dzm(nz)=dzm(nz-1);

% JGP NOTE: In MITgcm the rho grid lies exactly in the centers of the cells
% defined by the w grid but this is not the case for ROMS. So reuse HLS
% definition of midpoint depth.
                                                                   
dz_w            = -diff(z_w);                       % equivalent to dz(1:nz-1)=z(1:nz-1)-z(2:nz);  50 elements
zm_w            = z_w(1:end-1) + .5*dz_w;          % reuse HLS definition                         50 elements
dzm_w           = zeros(1,nz_w);                   % dzm is a weird object contructed for use     51 elements
dzm_w(2:nz_w-1) = zm_w(1:nz_w-2)-zm_w(2:nz_w-1);   %    in the ODE solver, so do it as per HLS
dzm_w(1)        = dzm_w(2);
dzm_w(nz_w)     = dzm_w(nz_w-1);
dzm_w = abs(dzm_w);

%        get dynamic modes
A = zeros(nz_w,nz_w);
B = zeros(nz_w,nz_w);
%             create matrices   
for i=2:nz_w-1
  A(i,i)   =  1/(dz_w(i-1)*dzm_w(i))  + 1/(dz_w(i)*dzm_w(i));
  A(i,i-1) = -1/(dz_w(i-1)*dzm_w(i));
  A(i,i+1) = -1/(dz_w(i)*dzm_w(i));
end
for i=1:nz_w
  B(i,i)=N2_w(i);
end
%             set boundary conditions
A(  1 ,1)    = -1.;
A(nz_w,1) = -1.;

[wmodes_w,e] = eig(A,B);

%          extract eigenvalues
e=diag(e);
%
ind=find(imag(e)==0);
e=e(ind);
wmodes_w=wmodes_w(:,ind);
%
ind=find(e>=1.e-10);
e=e(ind);
wmodes_w=wmodes_w(:,ind);
%
[e,ind]=sort(e);
wmodes_w=wmodes_w(:,ind);

nm=length(e);
ce=1./sqrt(e);

% JGP NOTE: At this point we have the w modes on the w grid. The top row
%   (z-0) is all zeros, as is the bottom row (z=H). Put the w modes onto
%   the rho grid

wmodes_r = nan(nz_r,nm);
pmodes_r = nan(nz_r,nm);

% JGP NOTE: the pmodes are renormalized in ROMS_calc_psi.m, so don't bother
% with it here.

for nn=1:nm
     wmodes_r(:,nn) = interp1(z_w, wmodes_w(:,nn), z_r, 'spline');
%      fig(1);clf;plot(-z_w,wmodes_w(:,nn));hold on;plot(-z_r,wmodes_r(:,nn),'r')
     pmodes_r(:,nn) = diff(wmodes_w(:,nn))./diff(z_w)';
%      fig(2);clf;plot(-z_r,pmodes_r(:,nn))
end;
wmodes_r(1,:)   = wmodes_w(2,:)/2;
wmodes_r(end,:) = wmodes_w(end-1,:)/2;

aaa=5;





% % % for i=1:nm
% % % %           calculate first deriv of vertical modes
% % %   pr=diff(wmodes(:,i));   
% % %   pr(1:nz-1)= pr(1:nz-1)./dz(1:nz-1)';
% % %   pr=pr*rho0*ce(i)*ce(i);
% % % %       linearly interpolate back to original depths
% % %   pmodes(2:nz-1,i)=.5*(pr(2:nz-1)+pr(1:nz-2));
% % %   pmodes(1,i)=pr(1);
% % %   pmodes(nz,i)=pr(nz-1);
% % % end




%keyboard
pout = [0 cumsum(dz_w)];

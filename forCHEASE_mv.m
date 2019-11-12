function [funnetcdf,fname_out,globalsvalues, namelist_struct]=forCHEASE_mv(fname, funnetcdf, t0, ichease, tension)
%
% gets equilibrium information at asked time
%
% takes output from netcdf output of TRANSP and write EXPEQ_name_time file needed for interface
% with CHEASE, containing the minimum information needed to reconstruct an equilibrium
%
% varargin{1}: directory name where netcdf file resides
% varargin{2}: file name
% varargin{3}: netcdf function
%      if varargin{1}, {2} or {3} are empty, ask for it
% varargin{4}: rr(time,points,rho)
% varargin{5}: yy(time,points,rho)
%      if rr, yy, not given, compute plasma boundary from funnetcdf and asymmetric moments
%

funnetcdf=funnetcdf.allvars;
timearr=funnetcdf.TIME.data;
pcur=funnetcdf.PCUR.data;
if isempty(t0)
    figure;grid on; zoom on
    plot(timearr,pcur)
    title('Ip vs time')

    t0=input('enter time for equilibrium interface: ');
    close(gcf)
end
[~, ii]=min(abs(timearr-t0));
itime=ii;
disp(['closest time chosen: ' num2str(timearr(itime))]);
t0=timearr(itime);

% plasma boundary
% if ~isempty(varargin) && length(varargin)>=4
%   % rr, yy given
%   R=varargin{4}(itime,:,end)./100;
%   Z=varargin{5}(itime,:,end)./100;
% else
%   % get R,Z at itime
[R,Z]=plboundTRANSP(funnetcdf,itime,1);
% end

% create EXPEQ file, need R0 and B0 to normalize
% R, Z in m
R=R/100;
Z=Z/100;
Rmax=max(R);
Rmin=min(R);
Zmax=max(Z);
Zmin=min(Z);
R0=0.5*(Rmax+Rmin);
Z0=0.5*(Zmax+Zmin);
aminor=0.5*(Rmax-Rmin);
bminor=0.5*(Zmax-Zmin);
kappa=bminor/aminor;

% R*B0 with R in cm
bzxr=funnetcdf.BZXR.data;
B0=bzxr(itime)/R0/100;

ii=findstr(fname,'.');
EXPEQend=[fname(1:ii-1) '_t' num2str(t0)];
ffname = sprintf('/tmp/%s/EXPEQ_',getenv('USER'), EXPEQend);
fid=fopen(ffname,'w');
% aspect ratio
fprintf(fid,'   %.12f\n',aminor/R0);
% Z0 normalised to R0
fprintf(fid,'   %.12f\n',Z0/R0);

% profiles
a_ptowb=funnetcdf.PTOWB.data; % (vs time,X ! X is in mid points of XB)
a_pplas=funnetcdf.PPLAS.data;
a_X=funnetcdf.X.data;
a_XB=funnetcdf.XB.data;
a_PLFLX=funnetcdf.PLFLX.data;
a_GFUN=funnetcdf.GFUN.data;
% CUR not really j.B: a_CUR=funnetcdf{'CUR'}(:);
% To get <j.B>/R0/<B.grad phi>: 
% <j.B> given by PLJB (*mu0 1e4 in A/cm^2?) on X mesh (zone centered mesh as opposed to zone boundary mesh XB)
% <B.grad phi>=G(psi) <1/R**2> with G(psi) from GFUN*BZXR/100 and <1/R**2> from GR2I on X mesh
% Thus <j.B>/R0/<B.grad phi> [MKSA] = mu0 1e4 (PLJB*1e4) / (R0/100) / (GFUN*BZXR/100) / (1e4*GR2I)
a_PLJB=funnetcdf.PLJB.data;
a_GR2I=funnetcdf.GR2I.data;

X=a_X(:,itime);
XB=a_XB(:,itime);
PLFLX=a_PLFLX(:,itime);
GFUN=a_GFUN(:,itime).*bzxr(itime)./100; % in T*m
%CUR=a_CUR(itime,:).*bzxr(itime)./100.*100^2; % in A/m2 ?
PLJB = 4*pi*1e+1*a_PLJB(:,itime);
GR2I = 1e4 * a_GR2I(:,itime);

% all on XB mesh as psi on XB mesh
% ptowb=interpos(13,X,a_ptowb(itime,:),XB); % could impose 0 derivative at r/a=0
ptowb=interpos(13,cat(1, 0, X),cat(1,a_ptowb(1,itime), a_ptowb(:,itime)),XB,0.0,[1 0],[0 0]); % imposes 0 derivative in center
% pplas=interpos(13,X,a_pplas(itime,:),XB);
pplas=interpos(13,cat(1, 0, X),cat(1, a_pplas(1,itime), a_pplas(:,itime)),XB,0.0,[1 0],[0 0]); % 0 derivative in center
iplot=1;
if iplot
  figure
  plot(X,a_ptowb(:, itime));
  hold on
  plot(XB,ptowb,'r--')
  plot(X,a_pplas(:, itime),'c')
  plot(XB,pplas,'g--')
  legend('ptowb on X','ptowb on XB with 0 der.','pplas on X','pplas on XB 0 der.')
end

% ptowb seems total pressure, so use that
pedge=ptowb(end);
% pedge normalised by: mu0/B0^2
mu0=4.e-7*pi;
fprintf(fid,'   %.12f\n',pedge*mu0/B0^2);

% plasma boundary (normalized to R0)
% seems R, Z have many turns, so take just first turn, once passed twice the same value of atan2
ii=find((atan2(R(1:end-1),Z(1:end-1))-atan2(R(1),Z(1))).*(atan2(R(2:end),Z(2:end))-atan2(R(1),Z(1)))<0);
if length(ii)>1
  fprintf(fid,'   %d\n',ii(2));
  fprintf(fid,'   %.12f     %.12f\n',[R(1:ii(2)) ;Z(1:ii(2))]/R0);
else
  fprintf(fid,'   %d\n',length(R));
  fprintf(fid,'   %.12f     %.12f\n',[R ;Z]/R0);
end

% pprime and GGprime profiles on sqrt(psi) mesh
% allow for use of jpar as well
%tension=input('enter return or smoothing value (tension), default=1e-4: ');
if isempty(tension)
  tension=1e-4;
end

sout= 0:1/41:1;
psiout=sout.^2 .* PLFLX(end); % in Webers so that derivative with psi has correct dimension
% Warning: use mesh with x coordinate between 0 and 1 so tension has same realtive effect
% Thus divide by x(end) and rescale derivative
xscale=PLFLX(end);
[p,pprime]=interpos(13,PLFLX./xscale,ptowb,psiout./xscale,tension);
pprime=pprime./xscale;
[G,Gprime]=interpos(13,PLFLX./xscale,GFUN,psiout./xscale,tension);
Gprime=Gprime./xscale;
GGprime=G.*Gprime;

jdotnorm=100.* PLJB ./ R0 ./ GFUN ./ GR2I;
[jpar]=interpos(13,PLFLX./xscale,jdotnorm,psiout./xscale,tension);

if iplot
  figure
  subplot(2,1,1)
  plot(PLFLX,ptowb)
  hold on
  plot(psiout,p,'r--')
  legend('p','p_{it}')
  subplot(2,1,2)
  plot(psiout,pprime,'r')
  hold on
  [~,pprime0]=interpos(13,PLFLX./xscale,ptowb,psiout./xscale,0.0);
  plot(psiout,pprime0./xscale,'k')
  legend('pprime','pprime spline')
  
  figure
  subplot(2,1,1)
  plot(PLFLX,GFUN)
  hold on
  plot(psiout,G,'r--')
  legend('G','G_{fit}')
  subplot(2,1,2)
  plot(psiout,GGprime,'r')
  hold on
  [G0,Gprime0]=interpos(13,PLFLX./xscale,GFUN,psiout./xscale,0.0);
  plot(psiout,G0.*Gprime0./xscale,'k')
  legend('GGprime','GGprime spline')
  
  figure
  plot(PLFLX,jdotnorm)
  hold on
  plot(psiout,jpar,'r')
  legend('<J.B>/R0/<B.GRAD PHI>','fit')
end

fprintf(fid,'   %d\n',length(sout));
fprintf(fid,'   1\n'); % means pprim and GGprime are given
% s-mesh
fprintf(fid,'   %.12f\n',sout);
% pprime normalised by mu0*R0^2/B0 from MKSA
fprintf(fid,'   %.12f\n',pprime*mu0.*R0.^2/B0);
% GGprime normalised by 1/B0 from MKSA
fprintf(fid,'   %.12f\n',GGprime/B0);
% jpar normalised by mu0 R0/B0 from MKSA
fprintf(fid,'   %.12f\n',jpar.*mu0.*R0/B0);

% part necessary for equilibrium reconstruction finished in EXPEQ
% now add some addtional information to check the equilibrium

fprintf(fid,'\n');
fprintf(fid,'   %.12f  R02 [m]\n',R0);
fprintf(fid,'   %.12f  B0 [T]\n',B0);
fprintf(fid,'   %.12g  CURRT, => Ip [A] = %.12g\n',pcur(itime)*mu0/R0/B0,pcur(itime));
fprintf(fid,'   %.12f  Z0 [m]\n',Z0);

RAXIS=funnetcdf.RAXIS.data; RAXIS=RAXIS(itime);
fprintf(fid,'   %.12g  RAXIS [m]\n',RAXIS/100);
ZAXIS=funnetcdf.YAXIS.data; ZAXIS=ZAXIS(itime);
fprintf(fid,'   %.12g  ZAXIS [m]\n',ZAXIS/100);
Q0=funnetcdf.Q0.data; Q0=Q0(itime);
fprintf(fid,'   %.12g q0\n',Q0);
Q=funnetcdf.Q.data; Q=Q(:,itime);
fprintf(fid,'   %.12g  qedge\n',Q(end));
q95=interpos(13,PLFLX/PLFLX(end),Q,[0.95 0.95]);
fprintf(fid,'   %.12g  q95\n',q95(1));
qmin=min(Q);
fprintf(fid,'   %.12g  qmin\n',qmin(1));
lio2=funnetcdf.LIO2.data; lio2=lio2(itime);
fprintf(fid,'   %.12g  li\n',2.*lio2);
BTEQ=funnetcdf.BTEQ.data; BTEQ=BTEQ(itime);
fprintf(fid,'   %.12g  EQUILIBRIUM BETA(TOROIDAL)\n',BTEQ);
BTDIA=funnetcdf.BTDIA.data; BTDIA=BTDIA(itime);
fprintf(fid,'   %.12g  DIAMAGNETIC BETA(TOROIDAL)\n',BTDIA);
BPEQ=funnetcdf.BPEQ.data; BPEQ=BPEQ(itime);
fprintf(fid,'   %.12g  EQUILIBRIUM BETA(POLOIDAL)\n',BPEQ);
BPDIA=funnetcdf.BPDIA.data; BPDIA=BPDIA(itime);
fprintf(fid,'   %.12g  DIAMAGNETIC BETA(POLOIDAL)\n',BPDIA);

fprintf(fid,'   %.12g  kappa\n',kappa);
fprintf(fid,'   %.12g  aminor\n',aminor);

fprintf(fid,'\n');
fprintf(fid,'from TRANSP file: %s%s\n',fname);
fprintf(fid,'time chosen: %f\n',t0);
fprintf(fid,'tension for profiles (smoothing): %g\n',tension);
fprintf(fid,'date: %s\n',date);
fclose(fid);

%ichease=0;
%ichease=input('Run chease? (0=no, 1=yes) ');

if ichease
  addpath('/home/osauter/chease/matlab','-end')
  [~,~,namelist_struct] = run_chease;
  namelist_struct.ncscal = 4;
  [fname_out,globalsvalues,namelist_struct,~] = run_chease(1, ffname);
end
return
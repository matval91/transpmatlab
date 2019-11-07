function [R, Z]=plboundTRANSP(funnetcdf,time_indx,bnd_only,varargin);
%
% calculates plasma boundary using asymmetric moments from TRANSP netcdf variables in funnetcdf
%
% if time is a single index, compute R,Z at this index
% if time is an array of indices, compute R(1:length(time_indx),...) and Z(1:length(time_indx),...)
%
% if bnd_only=1: compute only outside boundary => R(1:nbpoints) or R(time_indx,1:nbpoints)
% if bnd_only=0: compute all surfaces => R(1:nbpoints,flux_surface_XB) or ...
%                                        R(1:length(time_indx),1:nbpoints,flux_surface_XB)
%
% In this way, all the plasma booundaries are obtained with:
%    [R,Z]=plboundTRANSP(funnetcdf,1:length(time),1);
% or at given index only:
%    [R,Z]=plboundTRANSP(funnetcdf,14,1);
%

% get data from funnetcdf

% get Moments ASYMETRIC
for i=2:9
  c=int2str((i-1));
  mom(:,:,1,i)=funnetcdf.(strcat('RMC0',c)).data;
  mom(:,:,2,i)=funnetcdf.(strcat('RMS0',c)).data;
  mom(:,:,3,i)=funnetcdf.(strcat('YMC0',c)).data;
  mom(:,:,4,i)=funnetcdf.(strcat('YMS0',c)).data;
end
mom(:,:,1,1)=funnetcdf.('RMC00').data;
mom(:,:,3,1)=funnetcdf.('YMC00').data;
mom(:,:,2,1)=0;
mom(:,:,4,1)=0;

xx1=61;
inx=length(mom(:,1,1,1));
theta =linspace(0,2.*pi,xx1);
ms = 0:8;
co = cos(ms' * theta);
s =  sin(ms' * theta);

ix_eff=1:inx;
if bnd_only
  ix_eff=inx;
end

% in case time_indx not from 1:N
ii_time_indx(1:length(time_indx))=time_indx;
zrr=zeros(length(ii_time_indx),xx1,inx,9);
zyy=zeros(length(ii_time_indx),xx1,inx,9);
% for a few times (if length(time)<xx1)
if length(ii_time_indx)<xx1
  for itime_ii=1:length(ii_time_indx)
    itime=ii_time_indx(itime_ii);
    for im2=1:9
      for ix=ix_eff
         %keyboard;
        zrr(itime_ii,1:xx1,ix,im2)=mom(ix, itime,1,im2).*....
            co(im2,1:xx1)+mom(ix, itime,2,im2).*s(im2,1:xx1);
        zyy(itime_ii,1:xx1,ix,im2)=mom(ix, itime,3,im2).*...
            co(im2,1:xx1)+mom(ix, itime,4,im2).*s(im2,1:xx1);
      end
    end
  end
else
  for im2=1:9
    for ix=ix_eff
      for ith=1:xx1
        zrr(1:length(time_indx),ith,ix,im2)=mom(time_indx,ix,1,im2).*co(im2,ith) + ...
            mom(time_indx,ix,2,im2).*s(im2,ith);
        zyy(1:length(time_indx),ith,ix,im2)=mom(time_indx,ix,3,im2).*co(im2,ith) + ...
            mom(time_indx,ix,4,im2).*s(im2,ith);
      end
    end
  end
end
zrr2=sum(zrr,4);
zyy2=sum(zyy,4);

if length(ii_time_indx)==1
  if bnd_only
    R=zrr2(1,:,end);
    Z=zyy2(1,:,end);
  else
    R(:,:)=zrr2(1,:,:);
    Z(:,:)=zyy2(1,:,:);
  end
else
  if bnd_only
    R(ii_time_indx',1:xx1)=zrr2(:,:,end);
    Z(ii_time_indx',1:xx1)=zyy2(:,:,end);
  else
    R(ii_time_indx',1:xx1,1:inx)=zrr2(:,:,:);
    Z(ii_time_indx',1:xx1,1:inx)=zyy2(:,:,:);
  end
end  


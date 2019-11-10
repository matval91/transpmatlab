function nubeam_plot(transp, t)
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

%% collect dimensions
time = transp.coords.TIME.data;
rho = transp.coords.X.data; %time and index dependant
if isempty(t) || t<min(time) || t> max(time)
    fprintf('Set time to 0, timerange = %s - %s \n',num2str(min(time)), num2str(max(time)));
    t=0;
end
[~,ind]  = min(abs(t-time));

%% collect 1D data
p_inj = transp.allvars.PINJ.data*1e-6; % total injected power
p_ST = transp.allvars.BPSHI.data*1e-6; % shine-through
p_OL = transp.allvars.BPLIM.data*1e-6; % orbit losses
p_CX = transp.allvars.BPCXI.data*1e-6+transp.allvars.BPCXX.data*1e-6; % charge-exchange
p_e = transp.allvars.BPTE.data*1e-6;
p_i = transp.allvars.BPTI.data*1e-6;
p_th = transp.allvars.BPTH.data*1e-6;

neut = transp.allvars.NEUTT.data;
neut_DD = transp.allvars.NEUTX_DD.data;
neut_thnuc = transp.allvars.NEUTX.data;

%% collect 2D data
area = transp.allvars.DAREA.data; %differential area (area of one flux tube)
vol = transp.allvars.DVOL.data;%differential volume (volume of one flux tube)
j_beam = transp.allvars.CURB.data*10.; % kA/m^2
pe_beam = transp.allvars.PBE.data; %MW/m3
pi_beam = transp.allvars.PBI.data;%MW/m3
n_beam = transp.allvars.BDENS.data*1e6; %1/m3
p_beam = transp.allvars.PMHDF_IN.data;
%% integrate 2D data
I_beam = dot(j_beam/10.,area)*1e-3; %kA

%% plot 1D
figure('Position', [10 10 1000 1000]);
ax1=subplot(3,2,[1,3]);
hold on;
p=p_inj; plot(time, p, 'Color', 'k', 'DisplayName', 'Inj.');
p=p-p_ST; plot(time, p, 'Color', 'c', 'DisplayName', 'Shine-through');
p=p-p_OL; plot(time, p, 'Color', 'b', 'DisplayName', 'Orbit losses');
p=p-p_CX; plot(time, p, 'Color', 'g', 'DisplayName', 'CX');
p=p-p_e; plot(time, p, 'Color', 'b', 'DisplayName', 'electrons');
p=p-p_i; plot(time, p, 'Color', 'r', 'DisplayName', 'ions');
p=p-p_th; plot(time, p, 'Color', 'm', 'DisplayName', 'Thermalized');
hold off;
title(sprintf('TRANSP: %s',num2str(transp.shot)));
xlabel('t [sec]'); ylabel('P [MW]'); grid on; box on; legend show;

ax5=subplot(3,2,5);
hold on;
plot(time, neut,'k');
% plot(time, neut_DD,'r', 'DisplayName', 'DD neutrons');
% plot(time, neut_thnuc,'r', 'DisplayName', 'thermonuclear neutrons');
xlabel('t [sec]'); ylabel('neutrons [1/sec]'); grid on; box on; 
hold off;
linkaxes([ax1,ax5], 'x');

%% plot 2D

ax2=subplot(3,2,2);
hold on;
plot(rho(:,ind), pe_beam(:,ind), 'r', 'DisplayName', 'electrons')
plot(rho(:,ind), pi_beam(:,ind), 'k', 'DisplayName', 'ions')
hold off;
title(sprintf('time %s',num2str(time(ind))));
xlabel('\rho_{tor}'); ylabel('p [MW/m^3]'); grid on; box on; legend show;

ax4=subplot(3,2,4);
hold on;
plot(rho(:,ind), j_beam(:,ind), 'k')
hold off;
xlabel('\rho_{tor}'); ylabel('j [kA/m^2]'); grid on; box on;

ax6=subplot(3,2,6);
hold on;
plot(rho(:,ind), n_beam(:,ind), 'k')
hold off;
xlabel('\rho_{tor}'); ylabel('n [1/cm^3]'); grid on; box on;
linkaxes([ax2,ax4, ax6], 'x');

return

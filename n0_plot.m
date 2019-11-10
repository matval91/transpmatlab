function n0_plot(transp, t)
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);
%
% gets and plots informations about neutrals

%% collect dimensions
time = transp.coords.TIME.data;
rho = transp.coords.X.data; %time and index dependant
if isempty(t) || t<min(time) || t> max(time)
    fprintf('Set time to 0, timerange = %s - %s \n',num2str(min(time)), num2str(max(time)));
    t=0;
end
[~,ind]  = min(abs(t-time));

%% collect 2D data
% area = transp.allvars.DAREA.data; %differential area (area of one flux tube)
% vol = transp.allvars.DVOL.data;%differential volume (volume of one flux tube)
v_source = transp.allvars.DN0VD.data;
w_source = transp.allvars.DN0WD.data;

CXfastn = transp.allvars.N0BCXD0.data;
first_fastn = transp.allvars.N0BD0.data;
n0fast = first_fastn+CXfastn;
tot_source = v_source+w_source+n0fast;

W_temp = transp.allvars.T0WD.data;
V_temp = transp.allvars.T0VD.data;
%% integrate 2D data


%% plot 1D
figure('Position', [10 10 1000 1000]);
ax1=subplot(3,2,[1,3]);
hold on;
plot(time, tot_source(end,:)*1e6, 'k','DisplayName', 'Edge');
plot(time, tot_source(1,:)*1e6, 'r','DisplayName', 'Core');
hold off;
title(sprintf('TRANSP: %s',num2str(transp.shot)));
xlabel('t [sec]'); ylabel('n0 [1/m^3]'); grid on; box on; legend show;

ax5=subplot(3,2,5);
hold on;
%plot(time, neut,'k');
% plot(time, neut_DD,'r', 'DisplayName', 'DD neutrons');
% plot(time, neut_thnuc,'r', 'DisplayName', 'thermonuclear neutrons');
xlabel('t [sec]'); ylabel('neutrons [1/sec]'); grid on; box on; 
hold off;
linkaxes([ax1,ax5], 'x');

%% plot 2D
ax2=subplot(3,2,2);
hold on;
plot(rho(:,ind), tot_source(:,ind)*1e6, 'r', 'DisplayName', 'Total D')
plot(rho(:,ind), w_source(:,ind)*1e6, 'k', 'DisplayName', 'Wall D')
plot(rho(:,ind), v_source(:,ind)*1e6, 'b', 'DisplayName', 'Volume D')
hold off;
title(sprintf('time %s',num2str(time(ind))));
xlabel('\rho_{tor}'); ylabel('n_0^D [1/m^3]'); grid on; box on; legend show;
set(gca, 'YScale', 'log')

ax4=subplot(3,2,4);
hold on;
plot(rho(:,ind), W_temp(:,ind), 'k', 'DisplayName', 'Wall')
plot(rho(:,ind), V_temp(:,ind), 'r', 'DisplayName', 'Volume')
hold off;
xlabel('\rho_{tor}'); ylabel('T_0 [eV]'); grid on; box on; legend show;

ax6=subplot(3,2,6);
 hold on;
% plot(rho(:,ind), n_beam(:,ind), 'k')
hold off;
xlabel('\rho_{tor}'); ylabel('n [1/cm^3]'); grid on; box on;
linkaxes([ax2,ax4, ax6], 'x');

return

function nubeam_plot(nubeammat)
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);
%
% plots informations about nubeam

%% collect dimensions
time = nubeammat.time;
rho = nubeammat.rho; %time and index dependant
ind=nubeammat.ind;
%% collect 1D data
p_inj = nubeammat.d1.p_inj; % total injected power
p_ST = nubeammat.d1.p_ST; % shine-through
p_OL = nubeammat.d1.p_OL; % orbit losses
p_CX = nubeammat.d1.p_CX; % charge-exchange
p_e = nubeammat.d1.p_e;
p_i = nubeammat.d1.p_i;
p_th = nubeammat.d1.p_th;
neut = nubeammat.d1.neut;

%% collect 2D data
area = nubeammat.d2.area; %differential area (area of one flux tube)
j_beam = nubeammat.d2.j_beam; % kA/m^2
pe_beam = nubeammat.d2.pe_beam; %MW/m3
pi_beam = nubeammat.d2.pi_beam;%MW/m3
n_beam = nubeammat.d2.n_beam; %1/m3
pr_beam = nubeammat.d2.pr_beam;
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
title(sprintf('TRANSP: %s',num2str(nubeammat.id)));
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
plot(rho, pe_beam, 'r', 'DisplayName', 'electrons')
plot(rho, pi_beam, 'k', 'DisplayName', 'ions')
hold off;
title(sprintf('time %s',num2str(time(ind))));
xlabel('\rho_{tor}'); ylabel('p [MW/m^3]'); grid on; box on; legend show;

ax4=subplot(3,2,4);
hold on;
plot(rho, j_beam, 'k')
hold off;
xlabel('\rho_{tor}'); ylabel('j [kA/m^2]'); grid on; box on;

ax6=subplot(3,2,6);
hold on;
plot(rho, n_beam, 'k')
hold off;
xlabel('\rho_{tor}'); ylabel('n [1/cm^3]'); grid on; box on;
linkaxes([ax2,ax4, ax6], 'x');

return

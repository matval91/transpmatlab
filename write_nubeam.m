function write_nubeam(nubeammat)
ind=nubeammat.ind;
time = nubeammat.time(ind);
rho = nubeammat.rho; %time and index dependant
area = nubeammat.d2.area; %differential area (area of one flux tube)
j_beam = nubeammat.d2.j_beam; % kA/m^2
pe_beam = nubeammat.d2.pe_beam; %MW/m3
pi_beam = nubeammat.d2.pi_beam;%MW/m3
n_beam = nubeammat.d2.n_beam*1e-19; %1/m3
pr_beam = nubeammat.d2.pr_beam*1e-3;
output=[rho, area, j_beam, pe_beam, pi_beam, n_beam, pr_beam];

header1=sprintf('t= %d s \n',time);
header2=sprintf('rho_tor area(m2) j(kA/m2) pe(MW/m3) pi(MW/m3) n(10e19/m3) pr(kPa) \n');

% open a file for writing
dirrr=sprintf('/tmp/%s',getenv('USER'));
fname=sprintf('%s/NUBEAM%s_t%f.dat', dirrr, nubeammat.id, time);
fid = fopen(fname, 'w');
fprintf('\n\n NUBEAM written on %s \n\n\n', fname);
% print a title, followed by a blank line
fprintf(fid, header1);
fprintf(fid, header2);
% print values in column order
% two values appear on each row of the file
fprintf(fid, '%f  %f %f  %f %f  %f %f\n', output');
fclose(fid);

return
function test_transp
%% Function to test transp matlab routines

%% read defined file
fname='/home/vallar/Matlab/transp/63704V32.CDF';
transpmat=cdf2mat(fname);
time=0.45;

%% get, plot, export nubeam file
nubeam_struct=nubeam_get(transpmat,time);
nubeam_plot(nubeam_struct); close all;
nubeam_write(nubeam_struct);
disp('Removing file');
system( sprintf('rm -rf %s', nubeam_struct.fname) );

%test export of 
return
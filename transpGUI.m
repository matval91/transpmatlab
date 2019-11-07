set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16)
%
% some trials with netcdf, using 50725c01.cdf as test
%
% Installed matlab netcdf toolbox from web (http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html#GUIDE)
%
% then
%% read transp file

[fname,pname]=uigetfile('*.cdf','Open NetCDF File');
matcdf_per_name = cdf2mat(pname,fname);

%% getting variable names
varnames = fieldnames(matcdf_per_name.allvars);
for i=1:length(varnames)
  varnames_labels{i} = matcdf_per_name.allvars.(varnames{i}).label;
end
coordnames = fieldnames(matcdf_per_name.coords);
for i=1:length(coordnames)
  coordnames_labels{i} = matcdf_per_name.coords.(coordnames{i}).label;
end
%% GUI
matcdf.hfig_main=figure;
set(matcdf.hfig_main,'position',[50 80 900 750])
matcdf.hvarnames=uicontrol('style','listbox','string',varnames_labels,'position',[20 205 800 510]);
matcdf.hvarnames_text=uicontrol('style','text','string','variables (dims): long name [units]','position',[20 720 300 15]);

matcdf.hcoord=uicontrol('style','listbox','string',coordnames_labels,'position',[20 20 600 160]);
matcdf.hcoord_text=uicontrol('style','text','string','coordinates (size)','position',[20 185 100 15]);

matcdf.h(1)=uicontrol('style','pushbutton','String','Display variable','position',[650 20 200 20], ...
    'Callback','aa=matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')})');
matcdf.h(2)=uicontrol('style','pushbutton','String','Plot variable','position',[650 50 200 20], ...
    'Callback','figure(matcdf.hfig_plot);plot(matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).data)');

matcdf.h(3)=uicontrol('style','pushbutton','String','Plot transposed var.','position',[650 80 200 20], ...
    'Callback','figure(matcdf.hfig_plot);plot(matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).data'')');
cback = ['figure(matcdf.hfig_plot);plot(matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dim{1},' ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).data); ', ...
                    'xlabel([matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dimname{1} '' '' ', ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dimunits{1}]);', ...
                    'ylabel([matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).name '' '' ', ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).units]);'];
matcdf.h(4)=uicontrol('style','pushbutton','String','Plot vs dim 1','position',[650 110 200 20], ...
          'Callback',cback);
cback = ['figure(matcdf.hfig_plot);plot(matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dim{2},' ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).data); ', ...
                    'xlabel([matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dimname{2} '' '' ', ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).dimunits{2}]);', ...
                    'ylabel([matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).name '' '' ', ...
                    'matcdf_per_name.allvars.(varnames{get(matcdf.hvarnames,''Value'')}).units]);'];
matcdf.h(5)=uicontrol('style','pushbutton','String','Plot vs dim 2','position',[650 140 200 20], ...
          'Callback',cback);
matcdf.h(6)=uicontrol('style','pushbutton','String','Plot NUBEAM','position',[650 170 200 20], ...
          'Callback','nubeam_plot(matcdf_per_name)');
matcdf.hfig_plot=figure;
%%

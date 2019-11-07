%
% some trials with netcdf, using 50725c01.cdf as test
%
% Installed matlab netcdf toolbox from web (http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html#GUIDE)
%
% then

[fname,pname]=uigetfile('*.cdf','Open NetCDF File');
funnetcdf=netcdf.open([pname fname],'nowrite');

% all data of all vars
allinfo=ncinfo([pname fname]);
allvarids=netcdf.inqVarIDs(funnetcdf);
allvarnames={allinfo.Variables(:).Name};
if length(allvarnames) ~= length(allvarids)
  allinfo
  error('problem with Variables, may be several groups')
  return
end

[varnames_sorted,ind_sort_varnames]=sort(allvarnames);

% to find a variable:
% strmatch('GFUN',allvarnames,'exact')

% then data of var{1} of name allvarnames{1} (='TIME') can be obtained from:
% time_data=funnetcdf{'TIME'}(:);
% or
% time_data=funnetcdf{allvarnames{1}}(:);

% coordinates variables:
% >>coord(funnetcdf)
%    TIME
%    TIME3
%    X
%    XB
%    THETA
%    EION
%    KSID
%    RMJSYM
%    RMAJM
%    ILDEN
%    IVISB
%    INTNC

% construct relevant data strtucture

coordnames=strtrim({allinfo.Dimensions(:).Name});
[coordnames_sorted,ind_sort_coordnames]=sort(coordnames);
fields_variables_to_copy = {'Name', 'Dimensions', 'Size', 'Datatype'};
for i=1:length(coordnames_sorted)
  matcdf.coords(i).name = coordnames_sorted{i};
  matcdf.coords(i).index_allvarnames = strmatch(matcdf.coords(i).name,allvarnames,'exact');
  matcdf.coords(i).index_varnames_sorted = strmatch(matcdf.coords(i).name,varnames_sorted,'exact');
  matcdf.coords(i).varid = allvarids(matcdf.coords(i).index_allvarnames);
  matcdf.coords(i).data = netcdf.getVar(funnetcdf,matcdf.coords(i).varid);
  for ij=1:length(fields_variables_to_copy)
    matcdf.coords(i).(fields_variables_to_copy{ij}) = allinfo.Variables(matcdf.coords(i).index_allvarnames).(fields_variables_to_copy{ij});
  end
  matcdf.coords(i).units = strtrim(allinfo.Variables(matcdf.coords(i).index_allvarnames).Attributes(1).Value);
  matcdf.coords(i).long_name = strtrim(allinfo.Variables(matcdf.coords(i).index_allvarnames).Attributes(2).Value);
  matcdf.coords(i).label = [matcdf.coords(i).name ' ' num2str(matcdf.coords(i).Size) ': ' matcdf.coords(i).long_name];
  matcdf_per_name.coords.(matcdf.coords(i).name) = matcdf.coords(i);
end

for i=1:length(varnames_sorted)
  matcdf.vars(i).name = varnames_sorted{i};
  matcdf.vars(i).index_allvarnames = strmatch(matcdf.vars(i).name,allvarnames,'exact');
  matcdf.vars(i).varid = allvarids(matcdf.vars(i).index_allvarnames);
  matcdf.vars(i).data = netcdf.getVar(funnetcdf,matcdf.vars(i).varid);
  for ij=1:length(fields_variables_to_copy)
    matcdf.vars(i).(fields_variables_to_copy{ij}) = allinfo.Variables(matcdf.vars(i).index_allvarnames).(fields_variables_to_copy{ij});
  end
  matcdf.vars(i).units = strtrim(allinfo.Variables(matcdf.vars(i).index_allvarnames).Attributes(1).Value);
  matcdf.vars(i).long_name = strtrim(allinfo.Variables(matcdf.vars(i).index_allvarnames).Attributes(2).Value);
  matcdf.vars(i).label = matcdf.vars(i).name;
  for j=1:length(matcdf.vars(i).Dimensions)
    ij = strmatch(matcdf.vars(i).Dimensions(j).Name,coordnames_sorted,'exact');
    matcdf.vars(i).dim{j} = matcdf.coords(ij).data;
    matcdf.vars(i).dimunits{j} = matcdf.coords(ij).units;
    matcdf.vars(i).dimname{j} = matcdf.coords(ij).name;
    if j==1
      matcdf.vars(i).label = [matcdf.vars(i).label '(' matcdf.coords(ij).name];
    else
      matcdf.vars(i).label = [matcdf.vars(i).label ',' matcdf.coords(ij).name];
    end
    if j==length(matcdf.vars(i).Dimensions)
      matcdf.vars(i).label = [matcdf.vars(i).label ')'];
    end
  end
  matcdf.vars(i).label = [matcdf.vars(i).label ': ' matcdf.vars(i).long_name ' [' matcdf.vars(i).units ']'];
  matcdf_per_name.allvars.(matcdf.vars(i).name) = matcdf.vars(i);
end
top_attr_names = {allinfo.Attributes(:).Name};
ij = strmatch('shot',top_attr_names,'exact');
matcdf_per_name.shot = allinfo.Attributes(ij).Value;

%
% save 63704V32_matcdf.mat matcdf_per_name
clear matcdf
netcdf.close(funnetcdf);

% break
% to see the list:
%%

% load 63704V32_matcdf.mat
%%

varnames = fieldnames(matcdf_per_name.allvars);
for i=1:length(varnames)
  varnames_labels{i} = matcdf_per_name.allvars.(varnames{i}).label;
end
coordnames = fieldnames(matcdf_per_name.coords);
for i=1:length(coordnames)
  coordnames_labels{i} = matcdf_per_name.coords.(coordnames{i}).label;
end
%%
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

matcdf.hfig_plot=figure;
%%

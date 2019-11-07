function transpmat = cdf2mat(pfname)
%
% some trials with netcdf, using 50725c01.cdf as test
%
% Installed matlab netcdf toolbox from web (http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html#GUIDE)
%
% then
funnetcdf=netcdf.open(pfname,'nowrite');

% all data of all vars

allinfo=ncinfo(pfname);

allvarids=netcdf.inqVarIDs(funnetcdf);
allvarnames={allinfo.Variables(:).Name};
if length(allvarnames) ~= length(allvarids)
  allinfo;
  error('problem with Variables, may be several groups')
  return
end

[varnames_sorted,~]=sort(allvarnames);

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
[coordnames_sorted,~]=sort(coordnames);
fields_variables_to_copy = {'Name', 'Dimensions', 'Size', 'Datatype'};
for i=1:length(coordnames_sorted)
  matcdf.coords(i).name = coordnames_sorted{i};
  matcdf.coords(i).index_allvarnames = strmatch(matcdf.coords(i).name,allvarnames,'exact');
  matcdf.coords(i).index_varnames_sorted = strmatch(matcdf.coords(i).name,varnames_sorted,'exact');
  matcdf.coords(i).varid = allvarids(matcdf.coords(i).index_allvarnames);
  if strcmp(allinfo.Variables(matcdf.coords(i).index_allvarnames).Datatype,'single')
    matcdf.coords(i).data = netcdf.getVar(funnetcdf,matcdf.coords(i).varid,'double');
  else
    matcdf.coords(i).data = netcdf.getVar(funnetcdf,matcdf.coords(i).varid);
  end
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
  if strcmp(allinfo.Variables(matcdf.vars(i).index_allvarnames).Datatype,'single')
    matcdf.vars(i).data = netcdf.getVar(funnetcdf,matcdf.vars(i).varid,'double');
  else
    matcdf.vars(i).data = netcdf.getVar(funnetcdf,matcdf.vars(i).varid);
  end
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
transpmat = matcdf_per_name;
netcdf.close(funnetcdf);

clear matcdf
return 


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
allnames={allinfo.Variables(:).Name};
if length(allvars) ~= length(allvarids)
  allinfo
  error('problem with Variables, may be several groups')
  return
end

[names_sorted,ind_sort]=sort(allnames);
keyboard
% to find a variable:
% strmatch('GFUN',allnames,'exact')

% then data of var{1} of name allnames{1} (='TIME') can be obtained from:
% time_data=funnetcdf{'TIME'}(:);
% or
% time_data=funnetcdf{allnames{1}}(:);

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


% to see the list:
hfig=figure;
set(hfig,'position',[50 80 550 600])
hnames=uicontrol('style','listbox','string',names_sorted,'position',[20 250 200 300]);
hnames_text=uicontrol('style','text','string','variables','position',[20 555 80 15]);

coord_main=coord(funnetcdf);
% disp(' coordinates    size')
for i=1:length(coord_main)
  coord_names{i}=name(coord_main{i});
  coord_size{i}=num2str(size(coord_main{i}));
  coord_string{i}=[coord_names{i} '      (' coord_size{i} ')'];
  % disp(coord_string{i})
end

hcoord=uicontrol('style','listbox','string',coord_string,'position',[20 80 200 130]);
hcoord_text=uicontrol('style','text','string','coordinates (size)','position',[20 215 100 15]);

h(1)=uicontrol('style','pushbutton','String','Dim of variable','position',[250 20 250 20], ...
    'Callback','size(funnetcdf{allnames{ind_sort(get(hnames,''Value''))}})');
h(2)=uicontrol('style','pushbutton','String','Plot variable','position',[250 50 250 20], ...
    'Callback','figure(12);plot(funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(:))');

h(3)=uicontrol('style','pushbutton','String','Plot tranversed var.','position',[250 80 250 20], ...
    'Callback','figure(12);plot(funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(:)'')');

% $$$
% $$$ h(5)=uicontrol('style','pushbutton','String',['Plot var vs X (if Dim: ' num2str(size(funnetcdf{'X'})) ')'], ...
% $$$       'position',[250 140 250 20], ...
% $$$     'Callback','figure(12);plot(funnetcdf{''X''}(:)'',funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(:)'')');
% $$$
h(4)=uicontrol('style','pushbutton','String',['Plot var vs coord'], ...
      'position',[250 110 250 20], ...
    'Callback','figure(12);plot(funnetcdf{coord_names{get(hcoord,''Value'')}}(:)'',funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(:)'')');

% plot(funnetcdf{'XB'}(:)',funnetcdf{allnames{ind_sort(get(hnames,'Value'))}}(:)')

h(5)=uicontrol('style','pushbutton','String',['Choose time'], ...
    'position',[250 140 125 20],'Callback', ...
    ['if (size(str2num(get(h(6), ''string''))) < 2);' ...
      '[a it] = min(abs(funnetcdf{''TIME''}(:)-str2num(get(h(6), ''string''))));' ...
      'else;it=find(funnetcdf{''TIME''}(:) > min(str2num(get(h(6),''string'')))' ...
      '& funnetcdf{''TIME''}(:) < max(str2num(get(h(6),''string'')))); end;' ...
      'figure(12);plot(funnetcdf{coord_names{get(hcoord,''Value'')}}(it, :)'',' ...
      'funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(it,:)'')']);

h(6)=uicontrol('style','edit','String',[ '[' num2str(funnetcdf{'TIME'}(1)) ...
                    '  ' num2str(funnetcdf{'TIME'}(end)) ']' ], ...
          'position',[375 140 125 20]);
%  ,'Callback', ...
%    'figure(12);plot(funnetcdf{''TIME''}(:),funnetcdf{allnames{ind_sort(get(hnames,''Value''))}}(:))');

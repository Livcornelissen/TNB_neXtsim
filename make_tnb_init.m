close all
fpth='/home/jes/03_science/03_output/03_Meetings_Workshops/2025_SASIP/TNB_neXtsim/';

% fnme = 'init_new.nc';
% 
% fle = [fpth fnme];
% ncinfo(fle)
% msk = ncread(fle,'data/mask');
% 
% imagesc(msk)
% 
% nm = size(msk);
% 
% ns_dit = round(nm(2)/10);
% ew_dit = round(nm(1)/2);
% 
% 
% msk(1:ew_dit,round(nm(2)/2)+[0:ns_dit])=0;
% 
% 
% figure;imagesc(msk)
% 
% ncwrite(fle,'data/mask',msk)


%% create forcing file for polynya model
% 
%rogrd = '/home/jes/08_data/028/preparation/028_grd_RSSM38km24_V1.nc';
%template_grd = [fpth 'init_terranova.nc'];
initf = [fpth 'init_terranova.nc'];
frcfl = [fpth 'forcing_polynya.nc'];
if exist(frcfl,'file')
    delete(frcfl)
end
%copyfile(template_grd,nxgrd);

% create groups
ncid = netcdf.create(frcfl,"netcdf4");
strncid = netcdf.defGrp(ncid,"structure");
metncid = netcdf.defGrp(ncid,"metadata");

% create time group within metadata group
time_met_id = netcdf.defGrp(metncid,"time");

dum_dim = netcdf.defDim(time_met_id,'dim_dum',1);
tmp1=netcdf.defVar(time_met_id,'formatted','NC_STRING',dum_dim);

netcdf.putAtt(time_met_id,tmp1,'format',"%Y-%m-%dT%H:%M:%SZ")

tmp2=netcdf.defVar(time_met_id,'time','NC_INT64',dum_dim);
%netcdf.putAtt(time_met_id,tmp2,'units','seconds since 2023-04-01T00:00:00Z')
netcdf.putAtt(time_met_id,tmp2,'units','seconds since 1970-01-01T00:00:00Z')
dtancid = netcdf.defGrp(ncid,"data");

% read dimensions from init file and define them in forcing file
initid = netcdf.open(initf);
group_id = netcdf.inqNcid(initid, 'data');

[ndims, nvars, ngatts, unlimdimid] = netcdf.inq(group_id);

for n=0:ndims-1 % x and y
    [dim_nme,dim_len] = netcdf.inqDim(group_id,n);

    if strcmpi(dim_nme,'xdim') | strcmpi(dim_nme,'ydim')
        disp([char(dim_nme) ' ' num2str(dim_len)]);
        eval([dim_nme '=' num2str(dim_len) ';'])
    end
    %netcdf.defDim(dtancid,dim_nme,dim_len(i));
end
%% add time dimension
time = 12;

% create and populate dimensions in forcing file

dim_nme={'xdim','ydim','time'};
dim_len=[xdim, ydim, time];
for i=1:length(dim_nme)
    dimid(i) = netcdf.defDim(dtancid,char(dim_nme(i)),dim_len(i));
end

% create variables 
var_dim =   {'yx',       'yx'        't'  ,'tyx'  ,'tyx'  ,'tyx'  ,'tyx'  ,'tyx' ,'tyx'    ,'tyx' ,'tyx'};
var_nme={'longitude','latitude','time','dew2m','lw_in','sw_in','pair','tair','wind_speed','u','v'};
for i=1:length(var_nme)
    vdm = [contains(char(var_dim(i)),'x')*1 contains(char(var_dim(i)),'y')*2 contains(char(var_dim(i)),'t')*3];
    vdm(vdm==0)=[];
    var_id(i)=netcdf.defVar(dtancid,char(var_nme(i)),'NC_DOUBLE',[dimid(vdm)]);
end



st_lon=162;
en_lon=172;
str_lon=(en_lon-st_lon)/xdim;

st_lat=-76;
en_lat=-74;
str_lat=(en_lat-st_lat)/ydim;

lon = [st_lon:str_lon:en_lon-str_lon];
lat = [st_lat:str_lat:en_lat-str_lat];

[LON,LAT]=meshgrid(lon,lat);

k = find(contains(var_nme,'longitude')==1);
netcdf.putVar(dtancid, var_id(k), LON');

k = find(contains(var_nme,'latitude')==1);
netcdf.putVar(dtancid, var_id(k), LAT');

yf = [0:ydim-1]/max(ydim-1)*2;

wy_u = (-cos(yf*pi)+1)/2;
wx_u = exp(-[0:5/(xdim-1):5])';

w_u = wx_u * wy_u;

wx_v = [zv(xdim-round(xdim/2)-1) (exp([0:5/(round(xdim/2)):5])-1)/20]';
wx_v(wx_v>1)=1;

wy_v = ov(ydim); 
w_v = wx_v * wy_v;

u_max = 30;
v_max = 10;

u_wnd = u_max*w_u;
v_wnd = v_max*w_v;

%u_wnd=u_wnd';
%v_wnd=v_wnd';

u_wnd = repmat(u_wnd,[1 1 time]);
v_wnd = repmat(v_wnd,[1 1 time]);

u_wnd(:,:,4:6)=0;
u_wnd(:,:,10:12)=0;
v_wnd(:,:,4:6)=0;
v_wnd(:,:,10:12)=0;



% create 

k = find(matches(var_nme,'u')==1);
netcdf.putVar(dtancid, var_id(k), u_wnd);

k = find(matches(var_nme,'v')==1);
netcdf.putVar(dtancid, var_id(k), v_wnd);

wnd_s = sqrt(v_wnd.^2+u_wnd.^2);

k = find(matches(var_nme,'wind_speed')==1);
netcdf.putVar(dtancid, var_id(k), wnd_s);

var = wnd_s*0;

var(:)=0;
k = find(matches(var_nme,'dew2m')==1);
netcdf.putVar(dtancid, var_id(k), var);

var(:)=315.637;
k = find(matches(var_nme,'lw_in')==1);
netcdf.putVar(dtancid, var_id(k), var);

var(:)=0;
k = find(matches(var_nme,'sw_in')==1);
netcdf.putVar(dtancid, var_id(k), var);

var(:)=101325;
k = find(matches(var_nme,'pair')==1);
netcdf.putVar(dtancid, var_id(k), var);

var(:)=-30;
k = find(matches(var_nme,'tair')==1);
netcdf.putVar(dtancid, var_id(k), var);

 tme_var=1680321600+[0:11]*4*3600;
 k = find(matches(var_nme,'time')==1);
 netcdf.putVar(dtancid, var_id(k), tme_var);

%netcdf.putVar(dtancid,)

netcdf.close(ncid);


function [data, lons, lats] =  map_something(f, ops);
%function [data] =  map_something(f, ops);
% ops can be: 
%  "switch_latlon"
%  "fun" ( will apply function to column 9), 
%  "flatten" will get rid of depth dimension: if flatten='sum' or flatten='avg' adds up; if flatten is a function, will perform that operation (e.g. can set ops.flatten=@(x) sum(x,3) to get max value)
%  "field" will use a user input field instead of col. 9
%  "cloud" = 1 will assume that data has format: [lon lat dep Mmin Mmax rate XX] instead of [lon0 lon1 lat0 lat1 ...]
%  "clim" color axis

if nargin<2 ops=[]; end

if isfield(ops,'cloud')
  cloud=ops.cloud;
else 
  cloud=0;
end

if cloud
  lons=unique(f(:,1)');
  lats=unique(f(:,2)');
  deps=unique(f(:,3)');
  mags=unique(mean(f(:,4:5)'));
  f(:,9)=f(:,6);
else
  lons=unique(mean(f(:,1:2)'));
  lats=unique(mean(f(:,3:4)'));
  deps=unique(mean(f(:,5:6)'));
  mags=unique(mean(f(:,7:8)'));
end

nlon=length(lons);
nlat=length(lats);
ndep=length(deps);
nmags=length(mags);

if isfield(ops,'field')
   f(:,9)=ops.field;
end

if (nmags==1);
  if (nargin>1 && isfield(ops,'switch_latlon'))
    grid=reshape(f(:,9), nlon, nlat, ndep);
    grid=grid';
  else
    grid=reshape(f(:,9), nlat, nlon, ndep);
  end
else  
  if (nargin>1 && isfield(ops,'switch_latlon'))
    grid=reshape(f(:,9), nmags, nlon, nlat, ndep);
    grid=squeeze(sum(grid,1));
    grid=grid';
  else
    grid=reshape(f(:,9), nmags, nlat, nlon, ndep);
    grid=squeeze(sum(grid,1));
  end
end

if isfield(ops,'flatten')
   if strcmp(ops.flatten,'sum');
      ndep=1;
      grid=sum(grid,3);
  elseif strcmp(ops.flatten,'avg');
      ndep=1;
      grid=mean(grid,3);
  else
      ndep=1;
      grid=ops.flatten(grid);
  end
end


data=squeeze(grid);

if (nargin>1 && isfield(ops,'fun'))
  grid=ops.fun(grid);
end

nsp1=ceil(sqrt(ndep));
nsp2=ceil(ndep/nsp1);

figure
for n=1:ndep
 subplot(nsp1, nsp2, n);
 imagesc(lons, lats, grid(:,:,n))
 set(gca,'ydir','normal')
 %shading flat

 daspect_map;
 title([num2str(deps(n)) ' km depth']);
 if (nargin>1 && isfield(ops,'clim'))
   caxis(ops.clim);
 end
end


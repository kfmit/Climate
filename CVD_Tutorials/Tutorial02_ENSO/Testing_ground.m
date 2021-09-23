%% Test Script for Matlab Tut 02 ENSA

% Sccript to Read in Data
ncname = "sst.mon.mean_1982.01-2020.12.nc";

sstout = ncread(ncname,"sst");
time   = ncread(ncname,"time");
lat    = ncread(ncname,"lat");
lon    = ncread(ncname,"lon");

format long % type 'help format' into command line
t1 = time(1);

sst = permute(sstout,[3,1,2]); % The second argument indicates that new order of the original dimensions.
size(sst);
%%
for m = 1:12  %loop over month
    climSST(m,:,:) = mean(sst(m:12:end,:,:));
end
%whos climSST;

[ntime,nlon,nlat] = size(sst);% Get size of sst
nyrs = round(ntime/12);% Determine # of years
sst_mon   = reshape(sst,12,nyrs,nlon,nlat);% Reshape to separate month and year
climSST2  = squeeze(mean(sst_mon,2));% Take mean along year dimension
nanmax(abs(climSST2(:)-climSST(:)))% Check to make sure the mean absolute difference is 0 (between the two methods)

%% Test Ground 1: Subplots of Month SST
figure(6),clf
for m = 1:12
    subplot(3,4,m)
    % pcolor(...
    % caxis(...
    % colorbar ...
    % title( ...
    pcolor(lon,lat,squeeze(mean(sst(m:12:end,:,:)))');shading flat
    cb = colorbar;
    cb.Label.String = 'sst';
    title(['Map of Month ' num2str(m) ' Mean Sea Surface Temperature'])
    set(gca, 'Layer', 'top');
    xlabel('Longitude')
    ylabel('Latitude')
    caxis([0 35])
end

%% Testing Ground 2: SST Anomalies

sst(1,:,:);       % SST for January 1961
climSST(1,:,:);   % average SST for all Januarys

% SST anomaly for January 1961:
ssta_1961_1 = sst(1,:,:)-climSST(1,:,:);
% SST anomaly for January 1962:
ssta_1962_1 = sst(13,:,:)-climSST(1,:,:);

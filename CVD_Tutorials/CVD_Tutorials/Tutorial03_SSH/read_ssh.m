function [ssh,sshlon,sshlat,sshtime] = read_ssh(region)
%% This function loads sea level anomalies for given region
if region=='PA'
    sshlon = ncread('./SSH_Tropical_Eastern_Pacific2008_2012.nc','longitude');
    sshlat = ncread('./SSH_Tropical_Eastern_Pacific2008_2012.nc','latitude');

    % read all year of data
    sshtime = [];
    ssh = [];
    for year = [2008 2013]
        filename = dir(['./SSH_Tropical_Eastern_Pacific', num2str(year), '*.nc']);
        sshtime = [sshtime;double(ncread(filename.name,'time')+datenum(1950,01,01,00,00,00))];
        ssh = cat(3,ssh,ncread(filename.name,'adt'));
    end
    
elseif region=='NA'
    sshlon = ncread('./SSH_North_Atlantic2008_2010.nc','longitude');
    sshlat = ncread('./SSH_North_Atlantic2008_2010.nc','latitude');
    
    % read all year of data
    sshtime = [];
    ssh = [];
    for year = [2008 2011 2014 2017]
        filename = dir(['./SSH_North_Atlantic', num2str(year), '*.nc']);
        sshtime = [sshtime; double(ncread(filename.name,'time')+datenum(1950,01,01,00,00,00))];
        ssh = cat(3,ssh,ncread(filename.name,'adt'));
    end
else
    disp('This input is not valid, try ''NA'' or ''PA'' ')
end


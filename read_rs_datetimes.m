function [dtimes, dtimes_right, lon, lat] = read_rs_datetimes( path2sondes, assimilationTIME )
% read_rs_datetimes reads radiosonde launch times from available file names.
% It does not read observations here; only dates and launch coordinates needed
% to build the assimilation windows.

% for testing
% path2sondes = p2rs;

% Check all matching files in the "rs/<experiment>" folder.
list = dir([path2sondes, '*_sonde_*']);
dtimes = nan(size(list,1),2);
dtimes_right = nan(size(list,1),1); % Used for validation without repeated "dir" calls.
lat = dtimes_right; lon = dtimes_right;
for k = 1 : size(list,1)
    datfile  = [list(k).folder,'/',list(k).name];
    [~, fn, ~] = fileparts(datfile);
    if strcmp(fn(end-4), 'T')
        dayinfo  = datfile(end-18:end-4);
        dt_model = datenum(dayinfo, 'yyyy-mm-ddTHHMM');
        dt_right = dt_model;
    else
        % Only HH is available in the file name.
        dayinfo  = datfile(end-16:end-4);
        dt_model = datenum(dayinfo, 'yyyy-mm-ddTHH');
        dt_right = dt_model;
    end
    dtimes(k,1) = dt_model - datenum(0,0,0,0,round(assimilationTIME/2),0);
    dtimes(k,2) = dt_model + datenum(0,0,0,0,round(assimilationTIME/2),0);
    dtimes_right(k) = dt_right;
    % Coordinates.
    % Latitude: 38.776 Longitude: -9.126 -> first line of the data file.
    fileID = fopen(datfile,'r');
    cs = textscan(fileID,'%s %f %s %f',1);
    fclose(fileID);
    lat(k) = cs{2}; % Lat
    lon(k) = cs{4}; % Long
end

return

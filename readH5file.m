function [IDxyzCOORD, azimuth, elevation, siwv, dt, perc] = readH5file(fn, sites, site2exc, d1, d2, ll0, verbose)
% readH5file loads GNSS slant water vapor observations from an HDF5 file.
% It filters by station and time window, converts times to MATLAB datenum,
% projects station coordinates, and returns sorted observation arrays.
% testing
% ll0 = [x0, y0]; IDxyzCOORD = idsite;

disp('Expected time resolution is 30s')

sec = 30;
totalS = numel(d1:datenum(0,0,0,0,0,sec):d2);

x0 = ll0(1);
y0 = ll0(2);

tStart = tic;

nS = size(sites,1);
perc = zeros(nS,1);

if isempty(site2exc)
    site2exc = '----';
end

% HDF5 time reference.
tref = datenum(2000,1,1,0,0,0);

% Convert the requested interval to seconds since 2000-01-01.
d1_sec = (d1 - tref) * 86400;
d2_sec = (d2 - tref) * 86400;

% Cell arrays avoid large NaN-filled matrices.
ID_cell    = cell(nS,1);
az_cell    = cell(nS,1);
el_cell    = cell(nS,1);
siwv_cell  = cell(nS,1);
dt_cell    = cell(nS,1);
site2exclude = zeros(nS,1);

for k = 1:nS

    site = sites(k,:);

    if ismember(lower(site), lower(site2exc), 'rows')
        site2exclude(k) = 1;
        if verbose
            disp(['Site ',upper(site),'(',num2str(k,'%02d'),') | excluded (see site2exc in namelist) !'])
        end
        continue
    end

    g = ['/',site];

    % Read time first, still in seconds.
    dt_sec = h5read(fn, [g,'/dt']);

    idx = find(dt_sec >= d1_sec & dt_sec <= d2_sec);

    if isempty(idx)
        if verbose
            disp(['Site ',upper(site),'(',num2str(k,'%02d'),') | no data in range ', ...
                datestr(d1,'dd/mm/yy HH:MM:SS'), ' to ', datestr(d2,'dd/mm/yy HH:MM:SS')])
        end
        continue
    end

    i0 = idx(1);
    m  = numel(idx);

    % If the data are time-continuous, this is much faster.
    dt_sec_win = h5read(fn, [g,'/dt'], [i0,1], [m,1]);
    az     = h5read(fn, [g,'/azimuth'],     [i0,1], [m,1]);
    el     = h5read(fn, [g,'/elevation'],   [i0,1], [m,1]);
    pw     = h5read(fn, [g,'/spwv'],        [i0,1], [m,1]);

    % Coordinates are read only once per station.
    lon = h5read(fn, [g,'/lon']);
    lat = h5read(fn, [g,'/lat']);
    alt = h5read(fn, [g,'/alt']);

    [xx, yy] = gauss(x0, y0, lon, lat);

    % Convert only the selected time window.
    dt_datenum = tref + double(dt_sec_win) / 86400;

    % Store in cell arrays.
    ID_cell{k} = single([ ...
        repmat(k, m, 1), ...
        repmat(single(xx), m, 1), ...
        repmat(single(yy), m, 1), ...
        repmat(single(alt), m, 1)]);

    az_cell{k}   = single(az);
    el_cell{k}   = single(el);
    siwv_cell{k} = single(pw);
    dt_cell{k}   = dt_datenum;

    % Percentage of available epochs.
    ds = numel(unique(dt_datenum));
    perc(k) = (ds * 100) / totalS;

    if verbose

        disp(['Site ',upper(site),'(',num2str(k,'%02d'),') | ', ...
            datestr(d1,'dd/mm/yy HH:MM:SS'), ' to ', ...
            datestr(d2,'dd/mm/yy HH:MM:SS'), ' | ', ...
            num2str(perc(k),'%4.2f'),'% '])
    end
end

if verbose
    disp('Concatenating arrays')
end

% Remove empty cells.
valid = ~cellfun(@isempty, dt_cell);

IDxyzCOORD  = vertcat(ID_cell{valid});
azimuth     = vertcat(az_cell{valid});
elevation   = vertcat(el_cell{valid});
siwv        = vertcat(siwv_cell{valid});
dt          = vertcat(dt_cell{valid});

if verbose
    disp('Sorting arrays')
end

[dt, ind] = sort(dt);

azimuth      = double(azimuth(ind));
elevation    = double(elevation(ind));
siwv         = double(siwv(ind));
IDxyzCOORD   = double(IDxyzCOORD(ind,:));

tEnd = toc(tStart);

if verbose
    disp(['... done in ',num2str(tEnd, '%.2f'),' sec '])
end
end

function [site, xyz] = infoH5file(fn, p, verbose)
% infoH5file reads station metadata from a GNSS HDF5 file.
% It returns station IDs and geographic coordinates, with optional plotting
% and verbose summaries for quick inspection.

hi = h5info(fn);

siteList = {};
xyzList  = [];

for k = 1:numel(hi.Groups)

    gname = hi.Groups(k).Name;   % Example: '/hksc'.

    % Accept only groups of the form '/ssss'.
    if strlength(string(gname)) ~= 5
        continue
    end

    dnames = string({hi.Groups(k).Datasets.Name});

    % Confirm that lon, lat, and alt datasets exist.
    if ~all(ismember(["lon","lat","alt"], dnames))
        if verbose
            disp(['Skipping group ', gname, ' | missing lon/lat/alt'])
        end
        continue
    end

    s = gname(2:end);

    lon = h5read(fn, [gname,'/lon']);
    lat = h5read(fn, [gname,'/lat']);
    alt = h5read(fn, [gname,'/alt']);

    siteList{end+1,1} = s;
    xyzList(end+1,:) = [lon, lat, alt];

end

site = char(siteList);
xyz  = xyzList;

if verbose
    for k = 1:size(site,1)
        disp(['Site ',upper(site(k,:)),'(',num2str(k,'%02i'),') | lon = ', ...
            num2str(xyz(k,1),'%.3f'), ' deg | lat = ', ...
            num2str(xyz(k,2),'%.3f'), ' deg | hgt = ', ...
            num2str(xyz(k,3),'%.3f'), ' m'])
    end

    disp(['min(lon) = ',num2str(min(xyz(:,1)),'%.6f'), ...
        ', min(lat) = ',num2str(min(xyz(:,2)),'%.6f')])
end

if p
    figure
    hold on

    scatter(xyz(:,1), xyz(:,2), 45, xyz(:,3), 'filled')

    cb = colorbar;
    cb.Label.String = 'Height (m)';

    colormap(turbo)

    for k = 1:size(site,1)
        text(xyz(k,1), xyz(k,2)+0.005, ...
            [site(k,:), '(',num2str(k),')'], ...
            'FontSize', 8, ...
            'Color', 'k')
    end

    xlabel('Longitude (deg.)')
    ylabel('Latitude (deg.)')

    title(['min(lon) = ',num2str(min(xyz(:,1)),'%.6f'), ...
           ', min(lat) = ',num2str(min(xyz(:,2)),'%.6f')])

    axis equal
    grid on
    box on
end

end

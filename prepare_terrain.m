function [Xp, Yp, levels, xyzCart, demSRTM, site_num2del] = ...
    prepare_terrain(xy0, rot, numVx, numVy, resVx, resVy, resVz, sites, xyz, sIDs, perc, p, verbose)
% prepare_terrain builds the terrain-following tomographic grid.
% It projects station coordinates, checks domain membership, estimates the
% surface from station heights and optional SRTM data, and returns voxel levels.
% xy0   : lower-left grid reference and projection origin [x0, y0, lon0, lat0].
% rot   : grid rotation around the origin (x0, y0).
% numVx : number of voxels in x.
% numVy : number of voxels in y.
% resVx : voxel spacing in x, in meters.
% resVy : voxel spacing in y, in meters.
% resVz : vertical voxel spacing vector, in meters.
% sites, xyz : station IDs and coordinates read by infoH5file.m.
% sIDs  : station numeric identifiers.
% perc  : data availability percentage from readH5file.m for the selected window.
% p = 1 : plot the 2D and 3D grids.
% p = 2 : create publication-oriented tomographic-grid plots.

% for testing 
% xy0 = orig; xyz = xyzGEO; sIDs = est;

if verbose
    disp('Preparing the terrain (and 3d grid)... ')
end
tStart = tic;

method = 'linear'; % SRTM interpolation method; currently not used downstream.
% Number of stations used to influence voxel altitude; used to search for "maxd".
maxkd = 4;

% By default, all stations are used to prepare the DEM.
if verbose, disp('Using all sites to prepare terrain (default)'); end

% XYZ coordinates.
x0   = xy0(1); y0   = xy0(2); % Domain starting longitude/latitude at the lower-left corner.
lon0 = xy0(3); lat0 = xy0(4); % Projection-center longitude/latitude.
[x0_, y0_] = gauss(lon0, lat0, x0, y0);
nl = length(resVz);

% Voxel centers.
Xm = nan(numVy,numVx);
Ym = Xm;
for k = 0 : numVx - 1
    Xm(:,k+1) = x0_ + resVx*k;
end
for k = 0 : numVy - 1
    Ym(k+1,:) = y0_ + resVy*k;
end
Ym = flipud(Ym);

% Wall coordinates.
% X
Xi = Xm(1,1)-resVx/2;
Xf = Xm(1,end)+resVx/2;
Xe = linspace(Xi, Xf, size(Xm,2)+1);
Xp = nan(size(Xm,1)+1, size(Xm,2)+1);
for k = 1 : size(Xm,1)+1
    Xp(k, :) = Xe;
end
% Y
Yi = Ym(1,1)+resVy/2;
Yf = Ym(end,1)-resVy/2;
Ye = linspace(Yi, Yf, size(Ym,1)+1);
Yp = nan(size(Ym,1)+1, size(Ym,2)+1);
for k = 1 : size(Ym,2)+1
    Yp(:, k) = Ye';
end

% Rotation (anti-clockwise)
if rot ~= 0
    xCenter = Xp(end,1); % x pivot point
    yCenter = Yp(end,1); % y pivot point
    Xp = Xp - xCenter;
    Yp = Yp - yCenter;
    rXp = nan(size(Xp)); rYp = rXp;
    for r = 1 : size(Xp,1)
        for c = 1 : size(Xp,2)
            rXp(r,c) = Xp(r,c)*cosd(rot) - Yp(r,c)*sind(rot);
            rYp(r,c) = Xp(r,c)*sind(rot) + Yp(r,c)*cosd(rot);
        end
    end
    Xp = rXp + xCenter; Yp = rYp + yCenter;
else
    xCenter = Xp(end,1); % x pivot point
    yCenter = Yp(end,1); % y pivot point
end

% Geographic to Cartesian coordinates.
xyzCart = nan(size(xyz, 1),3);
[xyzCart(:,1), xyzCart(:,2)] = gauss(lon0, lat0, xyz(:,1), xyz(:,2));
xyzCart(:,3) = xyz(:,3);

% Check whether stations are outside the domain.
% Build the boundary polygon.
Xb = [Xp(1,:) Xp(:,end)' fliplr(Xp(end,:)) flipud(Xp(:,1))'];
Yb = [Yp(1,:) Yp(:,end)' fliplr(Yp(end,:)) flipud(Yp(:,1))'];
% Select stations inside the polygon.
mask = inpolygon(xyzCart(:,1), xyzCart(:,2), Xb, Yb);
site_num2del = zeros(length(mask),1);
sdel = [];
for k = 1 : length(mask)
    if (k == 1) && (mask(k) == 0)
        sdel = ['''',upper(sites(k,:)),''''];
        site_num2del(k) = sIDs(k);
    elseif (k > 1) && (mask(k) == 0)
        sdel = [sdel,';''',upper(sites(k,:)),''''];
        site_num2del(k) = sIDs(k);
    end
end
site_num2del(site_num2del==0) = [];
if ~isempty(sdel)
    disp(['Remove sites: site2exc = [',sdel,'] | out of domain | optional!!!']);
end

% cp for plots (without mask)
xscp = xyzCart(:,1);
yscp = xyzCart(:,2);
zscp = xyzCart(:,3);
sIDs_cp = sIDs;
ssss_cp = sites;

% Build the boundary polygon.
Xb = [Xp(1,:) Xp(:,end)' fliplr(Xp(end,:)) flipud(Xp(:,1))'];
Yb = [Yp(1,:) Yp(:,end)' fliplr(Yp(end,:)) flipud(Yp(:,1))'];
% Select stations inside the polygon.
mask = inpolygon(xyzCart(:,1), xyzCart(:,2), Xb, Yb);
xs = xyzCart(mask,1);
ys = xyzCart(mask,2);
zs = xyzCart(mask,3);
ssss = sites(mask,:);
sIDs = sIDs(mask,1);
perc = perc(mask,1);
num = length(xs);

% Distance from each voxel until the minimum number of sites ("maxkd") is included.
% All calculations here are horizontal.
xx = Xp(:); yy = Yp(:);
d = nan(size(xx,1),1);
for ii = 1 : size(xx,1)
    dists = sqrt((xx(ii) - xs).^2 + (yy(ii) - ys).^2);
    d_ = mink(dists,maxkd);
    d(ii) = d_(maxkd) + 10; % Add 10 m to ensure the farthest selected site is included.
end
maxd = reshape(d, size(Xp));

if ( 1 ) 

    % DEM in geographic coordinates.
    dlon = m_to_lon(resVx, mean(xyz(:,2))); dlon_ = m_to_lon(resVx/2, mean(xyz(:,2))); % Reach the wall.
    dlat = m_to_lat(resVy, mean(xyz(:,2))); dlat_ = m_to_lat(resVy/2, mean(xyz(:,2)));

    % DEM limits.
    lon = nan(size(Xp,2),1);
    lat = nan(size(Xp,1),1);
    for k = 1 : size(Xp,2), lon(k) = (x0-dlon_) + dlon*(k-1); end
    for k = 1 : size(Xp,1), lat(k) = (y0-dlat_) + dlat*(k-1); end
    [lon, lat] = meshgrid(lon, lat);
    lat = flipud(lat);
    latlim = [min(lat(:))-0.5, max(lat(:))+0.5]; lonlim = [min(lon(:))-0.5, max(lon(:))+0.5]; % DEM limits.

    % Always save in the "srtm/" directory; it is created if needed.
    % If the download fails, retrieve the file manually from:
    % https://terrain.ardupilot.org/SRTM3/North_America/
    if abs(latlim(1))<60 || abs(latlim(2))<60
        D = readhgt([latlim, lonlim], 'srtm3', 'outdir', [pwd,'/srtm'], 'interp');
        vlon = D.lon;
        vlat = D.lat;
        zdem = double(D.z);
        [vlon, vlat] = meshgrid(vlon, vlat);
        % Convert all coordinates to Cartesian.
        [x, y] = gauss(lon0, lat0, vlon, vlat);
        xdem = reshape(x, size(vlat)); % For plotting only.
        ydem = reshape(y, size(vlat));
        % pcolor(xdem, ydem, zdem); shading flat; colorbar
        % Output arrays.
        demSRTM(:,1) = xdem(:);
        demSRTM(:,2) = ydem(:);
        demSRTM(:,3) = zdem(:);
    else
        xdem = nan; ydem = nan; zdem = nan;
        demSRTM(:,1) = nan;
        demSRTM(:,2) = nan;
        demSRTM(:,3) = nan;
    end
    n = size(Xp,1)-1;
    m = size(Xp,2)-1;
end

% DEM interpolated from station heights.
F = scatteredInterpolant(xs, ys, zs, 'nearest', 'nearest');
fromSITE1 = F(Xp, Yp);
% figure; surf(Xp, Yp, fromSITE1); colorbar; hold on; scatter3(xs, ys, zs, 'filled')

%*******************************************************************
% Terrain following vertical grid mode
%*******************************************************************
% Inverse distance weighting using station elevations.
% Stations must provide a reasonable representation of local topography.
% Higher pexp values give less weight to more distant stations.
pexp = 6;
fromSITE2 = nan(size(Xp));
for jj = 1 : size(Xp,2)
    for ii = 1 : size(Xp,1)
        dist = sqrt((Xp(ii,jj) - xs).^2 + (Yp(ii,jj) - ys).^2);
        idx = dist <= maxd(ii,jj);
        w = 1./dist(idx).^pexp;
        if sum(isinf(w)) > 0
            fromSITE2(ii,jj) = ys(isinf(w));
        else
            fromSITE2(ii,jj) = sum(w.*zs(idx))/sum(w);
        end
    end
end
% figure; surf(Xp, Yp, fromSITE2); colorbar; hold on; scatter3(xs, ys, zs, 'filled')

% Select the interpolation with the smallest error.
% fromSITE1 -> from scatteredInterpolant
Ferr1 = scatteredInterpolant(Xp(:), Yp(:), fromSITE1(:), 'linear', 'none');
zs1 = Ferr1(xs, ys);
err1 = zs1 - zs;
vies1 = mean(zs1) - mean(zs);
% fromSITE2 -> from Inverse distance weighting
Ferr2 = scatteredInterpolant(Xp(:), Yp(:), fromSITE2(:), 'linear', 'none');
zs2 = Ferr2(xs, ys);
err2 = zs2 - zs;
vies2 = mean(zs2) - mean(zs);

e1 = sqrt(max(err1)^2+mean(err1)^2+vies1^2);
e2 = sqrt(max(err2)^2+mean(err2)^2+vies2^2);
disp('1:scatteredInterpolant; 2:Inverse distance weighting')
disp(['err1|bias1 = ',num2str(e1,'%.2f'),'|',num2str(vies1,'%.2f'),' [m]; err2|bias2 = ',num2str(e2,'%.2f'),'|',num2str(vies2,'%.2f'), ' [m]'])
if e1 < e2
    surf_ = fromSITE1 - vies1;
else
    surf_ = fromSITE2 - vies2;
end
Fveri = scatteredInterpolant(Xp(:), Yp(:), surf_(:), 'linear', 'none');
zsv = Fveri(xs, ys);
vies3 = mean(zsv) - mean(zs);
disp(['bias correction verification : ',num2str(vies3,'%.2f'),' [m]'])
% figure; surf(Xp, Yp, surf_); colorbar; hold on; scatter3(xs, ys, zs, 'filled')

new_zCart = nan(num,1);
for k = 1:num
    if isnan(xs(k)); continue; end
    found = false;
    for jj = 1:m
        for ii = 1:n
            % ----------------------------------------------------
            % 1. Define the cell.
            % ----------------------------------------------------
            Xcell = [Xp(ii,jj), Xp(ii+1,jj), Xp(ii+1,jj+1), Xp(ii,jj+1)];
            Ycell = [Yp(ii,jj), Yp(ii+1,jj), Yp(ii+1,jj+1), Yp(ii,jj+1)];
            % ----------------------------------------------------
            % 2. Check whether the station is inside the cell.
            % ----------------------------------------------------
            inside = inpolygon(xs(k), ys(k), Xcell, Ycell);
            if inside
                found = true;
                % ------------------------------------------------
                % 3. Compute the cell plane.
                % ------------------------------------------------
                Zcell = [surf_(ii,jj), surf_(ii+1,jj), surf_(ii+1,jj+1), surf_(ii,jj+1)];
                A = [Xcell(:), Ycell(:), ones(4,1)];
                coef = A \ Zcell(:);
                a = coef(1); b = coef(2); c = coef(3);
                inDEM = a*xs(k) + b*ys(k) + c;
                % ------------------------------------------------
                % 4. Correct the altitude.
                % ------------------------------------------------
                diff = zs(k) - (inDEM + 1);
                if diff < 0
                    % Station below the surface: correct it.
                    new_zCart(k) = zs(k) - diff;
                else
                    new_zCart(k) = zs(k);
                end
                % ------------------------------------------------
                % 5. Check the upper voxel boundary.
                % ------------------------------------------------
                if new_zCart(k) > inDEM + resVz(1)
                    error(['Station ', upper(ssss(k,:)), '(', num2str(sIDs(k),'%02i'),') is outside base voxel'])
                end
                % ------------------------------------------------
                % 6. Debug output.
                % ------------------------------------------------
                if verbose
                    fprintf('Station %s(%02i): DEM=%7.2f | Z=%7.2f | corr=%7.2f\n', ...
                        upper(ssss(k,:)), sIDs(k), inDEM, zs(k), new_zCart(k)-zs(k));
                end
                break
            end
        end
        if found; break; end
    end
    if ~found
        warning('Station %s(%02i) not in any cell', upper(ssss(k,:)), sIDs(k));
    end
end
% figure; surf(Xp, Yp, surf_); colorbar; hold on; scatter3(xs, ys, new_zCart, 'filled')
xyzCart = [];
zs = new_zCart;
xyzCart(:,1) = xs;
xyzCart(:,2) = ys;
xyzCart(:,3) = zs;

% Station positions, already extracted once.
xs_all = xs; ys_all = ys;
for jj = 1:m
    for ii = 1:n
        % --------------------------------------------------------
        % 1. Cell.
        % --------------------------------------------------------
        Xcell = [Xp(ii,jj), Xp(ii+1,jj), Xp(ii+1,jj+1), Xp(ii,jj+1)];
        Ycell = [Yp(ii,jj), Yp(ii+1,jj), Yp(ii+1,jj+1), Yp(ii,jj+1)];
        % --------------------------------------------------------
        % 2. Stations inside the cell.
        % --------------------------------------------------------
        inside = inpolygon(xs_all, ys_all, Xcell, Ycell);
        ind = find(inside);
        % --------------------------------------------------------
        % 3. More than one station.
        % --------------------------------------------------------
        if numel(ind) > 1
            if verbose
                fprintf('Sites in same voxel: ');
                for k = 1:numel(ind)
                    ii_s = ind(k);
                    fprintf('%s(%02i, perc=%.2f) ', upper(ssss(ii_s,:)), sIDs(ii_s), perc(ii_s));
                end
                fprintf('\n');
            end
        end
    end
end

if p == 1
    figure('position', [400, 400, 850, 400]);
    subplot(1,2,1); hold on
    if ~isnan(xdem(1))
        imagesc(xdem(1,:)./1000, ydem(:,1)./1000, zdem./1000); colorbar;
    end
    scatter(Xp(:)./1000, Yp(:)./1000, 10);
    for k = 1 : size(Yp,2)
        line([Xp(1,k)./1000, Xp(end,k)./1000], [Yp(1,k)./1000, Yp(end,k)./1000], 'color', 'red')
    end
    for k = 1 : size(Yp,1)
        line([Xp(k,1)./1000, Xp(k,end)./1000], [Yp(k,1)./1000, Yp(k,end)./1000], 'color', 'red')
    end
    xlabel('X (km)'); ylabel('Y (km)'); title('SRTM3')
    % Add the sites.
    scatter(xs./1000, ys./1000, 30, zs./1000, 'filled', 'MarkerEdgeColor', 'k', 'Marker', '^')
    for k = 1 : size(xs, 1)
        text(xs(k)./1000, ys(k)./1000+1, ssss(k, :), 'color', 'w')
    end
    xlim([min(Xp(:)./1000), max(Xp(:)./1000)]); ylim([min(Yp(:)./1000), max(Yp(:)./1000)]);
    axis equal
    
    subplot(1,2,2); hold on
    %imagesc(Xp(1,:), Yp(:,1), surf_); colorbar
    pcolor(Xp./1000, Yp./1000, surf_./1000); shading flat; colorbar
    scatter(Xp(:)./1000, Yp(:)./1000, 10);
    for k = 1 : size(Yp,2)
        line([Xp(1,k)./1000, Xp(end,k)./1000], [Yp(1,k)./1000, Yp(end,k)./1000], 'color', 'red')
    end
    for k = 1 : size(Yp,1)
        line([Xp(k,1)./1000, Xp(k,end)./1000], [Yp(k,1)./1000, Yp(k,end)./1000], 'color', 'red')
    end
    xlabel('X (km)'); ylabel('Y (km)'); title(['Interpolation method = ',method])
    % Add the sites.
    scatter(xscp./1000, yscp./1000, 50, zscp./1000, 'filled', 'MarkerEdgeColor', 'k', 'Marker', '^')
    for k = 1 : size(xscp, 1)
        text(xscp(k)./1000, yscp(k)./1000+1, [ssss_cp(k, :),'(',num2str(sIDs_cp(k)),')'])
    end
    % Rotation origin.
    scatter(xCenter./1000, yCenter./1000, 100, '^k', 'filled')
    axis equal
end

% Vertical levels.
levels = nan(size(surf_,1),size(surf_,2),nl+1);
alt = surf_;
levels(:,:,1) = alt; % Surface voxels start with the mean altitude inside each voxel.
sigma = linspace(0,0.8,nl+1);
for k = 2 : nl+1
    alt = alt + resVz(k-1);
    interval = (max(alt(:)) - min(alt(:)));
    a = min(alt(:)) + interval*sigma(k);
    b = max(alt(:)) - interval*sigma(k);
    if b > a
        levels(:,:,k) = rescale(alt, a, b);
    else
        % From this point onward, the levels are constant.
        levels(:,:,k) = rescale(alt, a, a);
    end
    %levels(3,5,2:end)-levels(3,5,1:end-1)
end

if p == 1
    figure; hold on
    for k = 1 : nl+1
        mesh(Xp./1000, Yp./1000, levels(:,:,k)./1000)
    end
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)')
    title('Hybrid terrain following grid')
    scatter3(xs./1000, ys./1000, zs./1000, 50, 'filled')
    view(-50,12)
    
    figure; hold on
    mesh(Xp./1000, Yp./1000, levels(:,:,1)./1000)
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)')
    title('Hybrid terrain following grid - first level')
    scatter3(xs./1000, ys./1000, zs./1000, 50, 'filled')
    for k = 1 : size(xs,1)
        text(xs(k)./1000, ys(k)./1000, zs(k)./1000, ssss(k,:))
    end
    zlim([min(min(levels(:,:,1)./1000)), max(max(levels(:,:,1)./1000))*2])
    view(-50,12)
end

tEnd = toc(tStart);
if verbose
    disp(['... done in ',num2str(tEnd, '%.3f'),' sec '])
end

end

function dlat = m_to_lat(dy, alat)
%
% dy   = latitude difference in meters
% dlat = latidute difference in degrees
% alat = average latitude between the two fixes
% Reference: American Practical Navigator, Vol II, 1975 Edition, p 5
% https://gitlab.ecosystem-modelling.pml.ac.uk/Gunduz/fvcom-toolbox/tree/cabd3e32d74b33ef93043d5b9174349156b98622/utilities

rlat = alat * pi/180;
m = 111132.09 * ones(size(rlat)) - ...
    566.05 * cos(2 * rlat) + 1.2 * cos(4 * rlat);
dlat = dy ./ m;
end

function dlon = m_to_lon(dx, alat)
%
% dlon = longitude difference in degrees
% dx   = longitude difference in meters
% alat = average latitude between the two fixes
% https://gitlab.ecosystem-modelling.pml.ac.uk/Gunduz/fvcom-toolbox/tree/cabd3e32d74b33ef93043d5b9174349156b98622/utilities

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dlon = dx ./ p;
end

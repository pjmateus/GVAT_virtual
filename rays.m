function [M, P, C] = rays(tri, vert, IdSites, coordSites, az, el, Xp, Yp, levels, p, verbose)
% rays computes ray exit coordinates at the outer tomographic-domain walls.
% It uses station antenna coordinates, azimuth/elevation angles, and the
% triangulated domain boundary to find each ray-domain intersection.
%         tri : delaunayTriangulation output, used only for plots.
%        vert : domain vertices from triDomain.m.
%     IdSites : numeric station IDs, used only for plots.
%  coordSites : [M, P, C] station antenna coordinates in meters.
%          az : azimuth vector in degrees.
%          el : elevation vector in degrees.
%  Xp, Yp, levels : wall coordinates from prepare_terrain.m, in meters.
%           p : if 1, plot each station.
% [M, P, C] contain ray exit coordinates at the external domain walls.
%

% testing
% tri=triD; vert=vertD; IdSites=idsite(ind,1); coordSites=idsite(ind,2:4); az=azimuth(ind); el=elevation(ind);
% tri=triD; vert=vertD; IdSites=virtualID; coordSites=sitesCOORD; az=AZv; el=ELv; 
% tri=triD; vert=vertD; IdSites=id_(ind1); coordSites=coord(ind1,:); az=taz; el=tel; % slants

if verbose
    disp('Preparing rays (obtaining the coordinates of the intersection with the outer walls of the domain)...')
end
tStart = tic;  

x_min = min(Xp(:));
x_max = max(Xp(:));
y_min = min(Yp(:));
y_max = max(Yp(:));
xdist = x_max-x_min;
ydist = y_max-y_min;
dh    = sqrt(xdist^2 + ydist^2) + 10; % Ensure dh always exits the domain.

% Determine each ray using the maximum horizontal distance in the domain.
Mf = coordSites(:,1) + dh .* sind(az);
Pf = coordSites(:,2) + dh .* cosd(az);
dinc = dh./cosd(el);
Cf = coordSites(:,3) + sqrt(dinc.^2 - dh.^2);
% Zenith rays produce Inf values.
idx = isinf(Cf);
Cf(idx) = max(levels(:)) + 10;

num_rays = size(coordSites,1); % all stations 
coord = nan(size(coordSites,1),3);
kk = 1;
for k = 1 : num_rays
    % Use "lineType" = "segment".
    orig = [coordSites(k,1), coordSites(k,2), coordSites(k,3)];
    dir  = [Mf(k), Pf(k), Cf(k)] - orig;
    [intersect, ~, ~, ~, xcoor] = TriangleRayIntersection(orig, dir, vert(:,1:3), vert(:,4:6), vert(:,7:9), 'lineType', 'segment');
    if sum(intersect==1) > 1 || sum(intersect) == 0
        error(['There are ',num2str(size(coord,1)-num_rays),' rays with two (or more) intersections, some station below the surface or some problem with the triangles!!!'])
    end
    if sum(intersect==1) == 1
        coord(kk,:) = xcoor(intersect==1,:);
        kk = kk + 1;
    end
    if el(k) == 90
        % Handle zenith rays separately.
        % Entering this block replaces the kk entry from the previous if.
        coord(kk-1,1) = coordSites(k,1);
        coord(kk-1,2) = coordSites(k,2);
        coord(kk-1,3) = Cf(k);
    end
end
ind = isnan(coord(:,1));
coord(ind,:) = [];
M = coord(:,1);
P = coord(:,2);
C = coord(:,3); % Referenced to the triangle bases in "vert".

if p
    K  = convexHull(tri);
    ns = length(unique(IdSites));
    ids= unique(IdSites);
    for k = 1 : ns
        idx = find(IdSites == ids(k));
        figure; hold on
        trisurf(K,tri.Points(:,1),tri.Points(:,2),tri.Points(:,3),'FaceAlpha', 0.1, 'EdgeColor', 'none', 'LineStyle', 'none')
        xx0 = coordSites(idx,1);
        yy0 = coordSites(idx,2);
        zz0 = coordSites(idx,3);
        %scatter3(xx0, yy0, zz0)
        xx1 = M(idx);
        yy1 = P(idx);
        zz1 = C(idx);
        plot3([xx0'; xx1'],[yy0'; yy1'],[zz0'; zz1'], '-k', 'LineWidth', 0.5);
        xlabel('X'); ylabel('Y'); zlabel('Z'); title(['SITE: ',num2str(ids(k))])
        axis equal
        view(-22,29);
    end
end

tEnd = toc(tStart);
if verbose
    disp(['... done in ',num2str(tEnd, '%.3f'),' sec '])
end

return

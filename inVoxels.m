function [voxelid, dist_, diff] = inVoxels(tri, vert, indVoxel, xyzI, xyzF)
% inVoxels determines which tomographic voxels are crossed by a ray.
% It uses the voxel triangulation from triVoxels.m and ray start/end
% coordinates to compute the distance traveled inside each intersected voxel.
% tri      : output from triVoxels.m, used only for plotting.
% vert     : vertices that reconstruct all triangles for TriangleRayIntersection.m.
% indVoxel : voxel index following ind2sub.m order; see
%            https://www.mathworks.com/help/matlab/ref/ind2sub.html.
%            The same information is stored in the last column of tri{:,3}: (r,c,v,ind).
% xyzI     : [M, P, C] initial ray coordinates.
% xyzF     : [M, P, C] final ray coordinates, determined by rays.m.
%
% voxelid : ID of each voxel, following the order imposed by triVoxels.m.
%   dist_ : distance in meters traveled by the ray inside each voxel.
%    diff : difference between the total ray length and the sum of per-voxel
%           ray lengths. It may be nonzero depending on ray definition,
%           point count, and whether the ray exits the domain.

% testing
% xyzI = xyzCi; xyzF = xyzCv(1,:);

% Plot switch for testing only; do not expose as an argument.
p = 0;
n = size(tri,1);

% For "lineType" = "ray" or "line".
%direcao = xyzF - xyzI;                      % Direction vector.
%direcao_unitaria = direcao / norm(direcao); % Normalized direction vector.
%dir = direcao_unitaria;

% For "lineType" = "segment".
orig = [xyzI(1), xyzI(2), xyzI(3)];
dir  = xyzF - xyzI;

% Call TriangleRayIntersection.m with all voxel triangles.
[intersect, t, ~, ~, xcoor] = TriangleRayIntersection(orig, dir, vert(:,1:3), vert(:,4:6), vert(:,7:9), 'lineType', 'segment');
num = intersect==1; 
if sum(num)
    voxels = indVoxel(num);
    coord  = xcoor(num,:);
    [~, ind] = sort(t(num));
    voxels = voxels(ind);
    coord = coord(ind,:);
else
    error('A ray does not intercept any voxels in the tomographic domain !!!')
    % A station was almost certainly not removed.
    %figure; hold on
    %for k = 1 : size(tri,1)
    %    tetramesh(tri{k,1},'FaceAlpha', 0.1, 'EdgeColor', 'k', 'LineStyle', '-');
    %end
    %plot3([xyzI(1),xyzF(1)],[xyzI(2),xyzF(2)],[xyzI(3),xyzF(3)],'-r','LineWidth',2.5)
end

% For testing.
if p
    figure; hold on
    uv = sort(unique(voxels));
    for k = 1 : length(uv)
        tetramesh(tri{uv(k),1},'FaceAlpha', 0.1, 'EdgeColor', 'k', 'LineStyle', '-'); % Or "none" for both.
    end
    for k = 1 : length(voxels)
        scatter3(coord(k, 1), coord(k, 2), coord(k, 3), 30, 'k', 'filled')
    end
    plot3(coord(:,1),coord(:,2),coord(:,3),'-k','LineWidth',5)
    plot3([xyzI(1),xyzF(1)],[xyzI(2),xyzF(2)],[xyzI(3),xyzF(3)],'-r','LineWidth',2.5)
end

% 1) Normal case:
%   one voxel, the initial voxel containing the station, is detected once;
%   the remaining voxels are detected twice, at entry and exit.
% 2) Problem case [1]:
%   when the ray follows the boundary between two voxels, the "border"
%   condition in TriangleRayIntersection.m does not resolve it, and one
%   voxel is detected three times.
% 3) Problem case [2]:
%   one voxel, not containing the station, is detected only once. This can
%   occasionally happen at the domain exit; the exact reason is still unknown.
% 4) Problem case [3]:
%   one voxel is detected four times. This can happen when the ray intersects
%   a triangle while being almost parallel to it.

u = unique(voxels);
m = length(u);
countVOXEL = nan(m,1);
for k = 1 : m
    countVOXEL(k) = sum(voxels == u(k));
end
ind1 = sum(countVOXEL == 1);
ind2 = sum(countVOXEL == 2); % The station may be below the voxel minimum altitude.
ind3 = sum(countVOXEL == 3);
ind4 = sum(countVOXEL == 4);

% After implementing domain rotation, one station can remain below the voxel base
% even when prepare_terrain.m places it 1 m above the surface.
% Remove that extra point in the first voxel here.
if ind2 == m
    % idx contains the coordinate index with the lowest elevation.
    % Remove this index because it is probably below the voxel base.
    idx = find(min(coord(:,3)));
    voxels(idx)  = [];
    coord(idx,:) = [];
    % Repeat the process.
    u = unique(voxels);
    m = length(u);
    countVOXEL = nan(m,1);
    for k = 1 : m
        countVOXEL(k) = sum(voxels == u(k));
    end
    ind1 = sum(countVOXEL == 1);
    %ind2 = sum(countVOXEL == 2); % The station may be below the voxel minimum altitude.
    ind3 = sum(countVOXEL == 3);    
    ind4 = sum(countVOXEL == 4);
end

normal = 0;
if ind1 == 1 && ind3 == 0 && ind4 == 0
    % 1) Normal case.
    normal = 1;
elseif ind1 == 1 && ind3 == 1 && ind4 == 0
    % 2) Problem case [1].
    v = u(countVOXEL == 3);
    iv = find(voxels == v);
    cv = coord(iv,:);
    d1 = norm(cv(1,:) - cv(2,:)); d2 = norm(cv(2,:) - cv(3,:)); d3 = norm(cv(1,:) - cv(3,:));
    if d1 > 0.1
        voxels(iv(end),:) = [];
        coord(iv(end),:) = [];
    elseif d2 > 0.1
        voxels(iv(1),:) = [];
        coord(iv(1),:) = [];
    elseif d3 > 0.1
        voxels(iv(2),:) = [];
        coord(iv(2),:) = [];
    end
    normal = 1;
elseif ind1 > 1
    % 3) Problem case [2], no solution.
    normal = 0;
elseif ind4 >= 1
    % 4) Problem case [3], no solution.
    normal = 0;
end

% If exactly one index has a single position, the first voxel enters the distance loop.
if normal == 1
    % Distance calculation.
    dist = nan(m,1);
    for k = 1 : m
        idx = find(voxels == u(k));
        if length(idx) == 1
            % Distance between the station point and the first intersection.
            % Only one index can contain a single element; all others must contain two.
            dist(k) = norm(coord(idx, :) - xyzI);
        elseif length(idx) == 2
            % Remaining voxels have two indices.
            dist(k) = norm( coord(idx(2), :)  -  coord(idx(1), :) );
        end
    end

    % New array with the size of tri{:,1}, allowing parfor use when this script
    % is called. Otherwise, this allocation is unnecessary.
    voxelid = nan(n, 1);
    voxelid(1:m) = u;

    dist_ = nan(n, 1);
    dist_(1:m) = dist;

    % Check the sum of per-voxel distances against the ray length from xyzI
    % to the coordinate where the ray leaves the domain.
    diff = abs(norm(coord(end, :) - xyzI) - sum(dist));

    % A large diff means the ray crosses an area without voxels, for example
    % low-elevation rays starting at Cabo da Roca and passing below Serra de Sintra.
else
    voxelid = nan(n, 1);
    voxelid(1:m) = nan(m,1);

    dist_ = nan(n, 1);
    dist_(1:m) = nan(m,1);
    diff = nan;
end

return

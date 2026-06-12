function [tri, vert] = triDomain(Xp, Yp, levels, p)
% triDomain builds a 3D triangulation of the outer tomographic domain.
% The resulting triangle list is formatted for TriangleRayIntersection.m and
% is used by rays.m to find domain-wall intersections.

% Single voxel representing the full domain.
% Constant "m" extends the boundary so rays are guaranteed to leave the domain
% in inVoxels.m and rays.m. It is later compensated with "wallRAYSexitAJUST".
m = 10; 
xg(1,1) = Xp(1,1) - m;
xg(2,1) = Xp(end,1) - m;
xg(3,1) = Xp(end,end) + m;
xg(4,1) = Xp(1,end) + m;
xg(5,1) = Xp(1,1) - m;
xg(6,1) = Xp(end,1) - m;
xg(7,1) = Xp(end,end) + m;
xg(8,1) = Xp(1,end) + m;
xg(9,1) = mean(xg(1:8));
% yy
yg(1,1) = Yp(1,1) + m;
yg(2,1) = Yp(end,1) - m;
yg(3,1) = Yp(end,end) - m;
yg(4,1) = Yp(1,end) + m;
yg(5,1) = Yp(1,1) + m;
yg(6,1) = Yp(end,1) - m;
yg(7,1) = Yp(end,end) - m;
yg(8,1) = Yp(1,end) + m;
yg(9,1) = mean(yg(1:8));
% zz
zg(1,1) = min(levels(:)) - m;
zg(2,1) = min(levels(:)) - m;
zg(3,1) = min(levels(:)) - m;
zg(4,1) = min(levels(:)) - m;
zg(5,1) = max(levels(:)) + m;
zg(6,1) = max(levels(:)) + m;
zg(7,1) = max(levels(:)) + m;
zg(8,1) = max(levels(:)) + m;
zg(9,1) = mean(zg(1:8));

X = [xg, yg, zg];
tri = delaunayTriangulation(X);
xyz = X(1:8,:);

% The polygon has 12 triangles, two per face.
t = 12; % Do not change.
vert = nan(t, 9);
i1 = 1;

% Two lower triangles.
% 1
T1 = [ xyz(1,1), xyz(1,2), xyz(1,3) ];
T2 = [ xyz(2,1), xyz(2,2), xyz(2,3) ];
T3 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];
% 2
T4 = [ xyz(1,1), xyz(1,2), xyz(1,3) ];
T5 = [ xyz(4,1), xyz(4,2), xyz(4,3) ];
T6 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];

% Two upper triangles.
% 1
T7 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
T8 = [ xyz(6,1), xyz(6,2), xyz(6,3) ];
T9 = [ xyz(7,1), xyz(7,2), xyz(7,3) ];
% 2
T10 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
T11 = [ xyz(8,1), xyz(8,2), xyz(8,3) ];
T12 = [ xyz(7,1), xyz(7,2), xyz(7,3) ];

% Lateral triangles.
% Lower-left triangle.
T13 = [ xyz(1,1), xyz(1,2), xyz(1,3) ];
T14 = [ xyz(2,1), xyz(2,2), xyz(2,3) ];
T15 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
% Upper-left triangle.
T16 = [ xyz(2,1), xyz(2,2), xyz(2,3) ];
T17 = [ xyz(6,1), xyz(6,2), xyz(6,3) ];
T18 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
% Lower-right triangle.
T19 = [ xyz(4,1), xyz(4,2), xyz(4,3) ];
T20 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];
T21 = [ xyz(8,1), xyz(8,2), xyz(8,3) ];
% Upper-right triangle.
T22 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];
T23 = [ xyz(7,1), xyz(7,2), xyz(7,3) ];
T24 = [ xyz(8,1), xyz(8,2), xyz(8,3) ];
% Lower-north triangle.
T25 = [ xyz(1,1), xyz(1,2), xyz(1,3) ];
T26 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
T27 = [ xyz(4,1), xyz(4,2), xyz(4,3) ];
% Upper-north triangle.
T28 = [ xyz(5,1), xyz(5,2), xyz(5,3) ];
T29 = [ xyz(8,1), xyz(8,2), xyz(8,3) ];
T30 = [ xyz(4,1), xyz(4,2), xyz(4,3) ];
% Lower-south triangle.
T31 = [ xyz(2,1), xyz(2,2), xyz(2,3) ];
T32 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];
T33 = [ xyz(6,1), xyz(6,2), xyz(6,3) ];
% Upper-south triangle.
T34 = [ xyz(6,1), xyz(6,2), xyz(6,3) ];
T35 = [ xyz(7,1), xyz(7,2), xyz(7,3) ];
T36 = [ xyz(3,1), xyz(3,2), xyz(3,3) ];

% Format read by TriangleRayIntersection.m; keep all triangles together to avoid one call per triangle.
vert(i1, 1:3) = T1;
vert(i1, 4:6) = T2;
vert(i1, 7:9) = T3;
vert(i1+1, 1:3) = T4;
vert(i1+1, 4:6) = T5;
vert(i1+1, 7:9) = T6;
vert(i1+2, 1:3) = T7;
vert(i1+2, 4:6) = T8;
vert(i1+2, 7:9) = T9;

vert(i1+3, 1:3) = T10;
vert(i1+3, 4:6) = T11;
vert(i1+3, 7:9) = T12;
vert(i1+4, 1:3) = T13;
vert(i1+4, 4:6) = T14;
vert(i1+4, 7:9) = T15;
vert(i1+5, 1:3) = T16;
vert(i1+5, 4:6) = T17;
vert(i1+5, 7:9) = T18;

vert(i1+6, 1:3) = T19;
vert(i1+6, 4:6) = T20;
vert(i1+6, 7:9) = T21;
vert(i1+7, 1:3) = T22;
vert(i1+7, 4:6) = T23;
vert(i1+7, 7:9) = T24;
vert(i1+8, 1:3) = T25;
vert(i1+8, 4:6) = T26;
vert(i1+8, 7:9) = T27;

vert(i1+9, 1:3) = T28;
vert(i1+9, 4:6) = T29;
vert(i1+9, 7:9) = T30;
vert(i1+10, 1:3) = T31;
vert(i1+10, 4:6) = T32;
vert(i1+10, 7:9) = T33;
vert(i1+11, 1:3) = T34;
vert(i1+11, 4:6) = T35;
vert(i1+11, 7:9) = T36;

if p
    K = convexHull(tri);
    figure; hold on
    trisurf(K,tri.Points(:,1),tri.Points(:,2),tri.Points(:,3), 'FaceAlpha',0.3, 'EdgeColor', 'none', 'LineStyle', 'none')
    xlabel('X'); ylabel('Y'); zlabel('Z')
    axis equal
    view(-22,29);
end

return

function [tri, vert, indVoxel] = triVoxels(Xp, Yp, levels, p) 
% triVoxels triangulates every tomographic voxel.
% It reads the eight vertices of each voxel, preserves their required order,
% and builds triangle lists compatible with TriangleRayIntersection.m.
%
% Lower voxel face.
%           1            4
% Order: (i,j,k)   - (i,j+1,k)
%           |            |
%        (i+1,j,k) - (i+1,j+1,k)
%           2            3
%
% Upper voxel face.
%           5            8
% Order: (i,j,k)   - (i,j+1,k)
%           |            |
%        (i+1,j,k) - (i+1,j+1,k)
%           6            7
%
% PLOTTING ONLY: delaunayTriangulation builds a 3D triangulation for each tomographic voxel.
% Each 3D voxel triangulation is saved in a cell array.
% INPUT:
% Xp, Yp, levels : tomographic wall coordinates.
%               p : if 1, create plots, which is slower.
% OUTPUT:
% tri      : delaunayTriangulation cell array, used only to plot individual voxels.
%            For inVoxels.m testing, inspect this variable inside the loop.
% vert     : triangle vector processed by TriangleRayIntersection.m.
%            It contains all triangles for every voxel in the domain.
% indVoxel : voxel-index vector for use with ind2sub.m.
% 
if p, figure; hold on; end

% mem
xg = nan(9,1);
yg = xg;
zg = xg;

n = size(Xp, 1)-1; % Rows.
m = size(Xp, 2)-1; % Columns.
l = size(levels, 3)-1;
tri = cell(n*m*l, 3);

s = 1;
for k = 1 : l
    for j = 1 : m
        for i = 1 : n
            % xx    
            xg(1) = Xp(i,j);
            xg(2) = Xp(i+1,j);
            xg(3) = Xp(i+1,j+1);
            xg(4) = Xp(i,j+1);       
            xg(5) = Xp(i,j);
            xg(6) = Xp(i+1,j);
            xg(7) = Xp(i+1,j+1);
            xg(8) = Xp(i,j+1);   
            xg(9) = mean(xg(1:8));
            % yy       
            yg(1) = Yp(i,j);
            yg(2) = Yp(i+1,j);
            yg(3) = Yp(i+1,j+1);
            yg(4) = Yp(i,j+1);       
            yg(5) = Yp(i,j);
            yg(6) = Yp(i+1,j);
            yg(7) = Yp(i+1,j+1);
            yg(8) = Yp(i,j+1);     
            yg(9) = mean(yg(1:8));
            % zz   
            zg(1) = levels(i,j,k);
            zg(2) = levels(i+1,j,k);
            zg(3) = levels(i+1,j+1,k);
            zg(4) = levels(i,j+1,k);       
            zg(5) = levels(i,j,k+1);
            zg(6) = levels(i+1,j,k+1);
            zg(7) = levels(i+1,j+1,k+1);
            zg(8) = levels(i,j+1,k+1);      
            zg(9) = mean(zg(1:8));

            X = [xg, yg, zg];
            tri{s,1} = delaunayTriangulation(X); % Contains all voxel triangles, used in plots.
            tri{s,2} = X(1:8,:);                 % Contains all voxel vertices.
                                                 % Order: (i,j,k)   - (i,j+1,k)
                                                 %           |            |
                                                 %        (i+1,j,k) - (i+1,j+1,k)
                                                 % Plus level k+1.
            tri{s,3} = [i,j,k,s];                % Contains (row, col, vertical level, voxel index).
                                                 % Index "s" can be used with ind2sub.m; see https://www.mathworks.com/help/matlab/ref/ind2sub.html.
            if p
                K = convexHull(tri{s,1});
                trisurf(K,tri{s,1}.Points(:,1),tri{s,1}.Points(:,2),tri{s,1}.Points(:,3), 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'LineStyle', '-')
                xlabel('X'); ylabel('Y'); zlabel('Z');
                axis equal
                view(-22,29);
            end 
            s = s + 1;
        end
    end
end

% Each polygon has 12 triangles, two per face.
t = 12; % Do not change.
n = size(tri,1);
vert = nan(n*t, 9);
indVoxel = nan(n*t, 1);
i1 = 1;

for k = 1 : n
    xyz  = tri{k,2};

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

    indVoxel(i1:i1+t-1)  = k;
    i1 = i1 + t;
end

return

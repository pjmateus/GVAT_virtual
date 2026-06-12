function [AZ_, EL_, siteID, XYZsite, methodID] = virtual_rays(mpc, ids, Xp, Yp, levels, minEL, target, p, verbose)
% virtual_rays generates synthetic ray directions for GNSS tomography.
% It builds azimuth/elevation pairs using one or more targeting strategies,
% filters rays below the minimum elevation, and returns station/method metadata.
% mpc    : [M, P, C] station coordinates in meters.
% ids    : numeric IDs of active stations.
% Xp, Yp, levels : 3D tomographic grid.
% minEL  : minimum elevation; lower virtual rays are discarded relative to station height.
% target : vector enabling/disabling each method:
%          [voxelsWALL, fixAZandELcycle, centerOFallVOXELS, riseRAYSinSOMEvoxels, centerOFVOXELSrayMINnum, centerOFVOXELSrayMINnumV2]
%      p : 1-by-6 array enabling/disabling plots by method.

% FIRST LEVEL ONLY.
% At the first level, use the minimum altitude and subtract "drop" meters.
% This attempts to maximize first-level rays above the minimum elevation. If
% drop = 0, rays point to the voxel center. This may introduce discontinuities
% between the first two levels, so check small domains carefully.
drop = 10; % m

% for testing
% method = [0,0,0,1,0,0]; p = [0,0,0,1,0,0];
% mpc = xyzCart; ids = est; target = method; 

% If plotting is active for several methods, assign a different ray color to each.
% Do not use red because it marks removed rays.
colors = [0.9290, 0.6940, 0.1250;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880;
          0.3010, 0.7450, 0.9330;
          0.0000, 0.4470, 0.7410;   
          0.0000, 0.0000, 0.0000]; % p(5) and p(6) always create a new plot.

if sum(target) == 0
    error('You need to define at least one target, order: 1:voxelsWALL, 2:fixAZandELcycle, 3:centerOFallVOXELS, 4:riseRAYSinSOMEvoxels, 5:centerOFVOXELSrayMINnum, 6:centerOFVOXELSrayMINnumV2], use 0 or 1 to activate.')
end

% 1) 
% Point only to wall voxels.
voxelsWALL = false; 

% 2)
% Independent of voxels; use fixed AZ and EL values in a cycle.
fixAZandELcycle = false; 
azstep = 3; % deg
elstep = 3; % deg

% 3) 
% Point to the center of every voxel.
centerOFallVOXELS = false; 

% 4) 
% Increase the number of rays only at selected levels, while still targeting all voxels.
riseRAYSinSOMEvoxels= false; 
% If the levels are arrays, rays are increased in two domain slabs.
% Example: minLEVEL = [1, 12] and maxLEVEL = [3, 15] targets levels 1-3
% and 12-15. The upper slab also increases rays across the whole domain.
minLEVEL = 5; %[1, 12];
maxLEVEL = 10; %[3, 14];
numRAYS  = 4; % Number of rays per voxel: (numRAYS-1)^2.
                             
% 5)             
% Each voxel must have at least x rays from the x nearest stations.
% Note that each ray continues through other voxels.
centerOFVOXELSrayMINnum = false;   
                                   
% 6)               
% Improved version of the previous method.
% This method attempts to remove unnecessary rays.
centerOFVOXELSrayMINnumV2 = false;
maxNUMrays = 4; % Valid for centerOFVOXELSrayMINnum and centerOFVOXELSrayMINnumV2.

% Add zenith rays to the virtual set.
addZENITH = 1;

% Memory allocation and grid dimensions.
n = size(Xp, 1)-1; % Rows.
m = size(Xp, 2)-1; % Columns.
l = size(levels, 3)-1;
numberOFsites = length(ids);

% Terminal information.
% Options.
if target(1)
    voxelsWALL = true;
    if verbose, disp('voxelsWALL= 1'); end
else
    if verbose, disp('voxelsWALL= 0'); end
end

if target(2)
    fixAZandELcycle = true;
    stepAZ = 0 : azstep : 359;    
    stepEL = minEL : elstep : 89; 
    if verbose, disp(['fixAZandELcycle= 1 | stepAZ= [',num2str(stepAZ(1)),',',num2str(stepAZ(2)-stepAZ(1)),',',num2str(stepAZ(end)),'] | stepEL= [',...
            num2str(stepEL(1)),',',num2str(stepEL(2)-stepEL(1)),',',num2str(stepEL(end)),']']); end
else
    if verbose, disp('fixAZandELcycle= 0'); end
end

if target(3) == 1
    centerOFallVOXELS = true;
    if verbose, disp('centerOFallVOXELS= 1'); end
else
    if verbose, disp('centerOFallVOXELS= 0'); end
end

if target(4)
    riseRAYSinSOMEvoxels = true;
    if verbose
        if length(minLEVEL) == 1
            disp(['riseRAYSinSOMEvoxels = 1 | minLEVEL= ',num2str(minLEVEL),' | maxLEVEL= ',num2str(maxLEVEL), ' | numRAYS= ', num2str(numRAYS), ' (by voxel)']);
        else
            disp(['riseRAYSinSOMEvoxels = 1 | minLEVEL= [',num2str(minLEVEL(1)),',',num2str(minLEVEL(2)),...
                '] | maxLEVEL= [',num2str(maxLEVEL(1)),',',num2str(maxLEVEL(2)),'] | numRAYS= ', num2str(numRAYS), ' (by voxel)']);
        end
    end
end

if target(5)
    centerOFVOXELSrayMINnum = true;
    if verbose, disp('centerOFVOXELSrayMINnum= 1'); end
    % Minimum rays per voxel, equivalent to the nearest stations for that voxel.
    % A voxel must be visited "maxNUMrays" times, matching "maxNUMrays" stations.
    % Memory allocation.
    siteIDS = cell(n*m*l,1);
    sz = [n,m,l];
    for kk_ = 1 : l
        for ii_ = 1 : m
            for jj_ = 1 : n
                % Voxel center.
                xxC = (Xp(jj_,ii_) + Xp(jj_+1,ii_) + Xp(jj_+1,ii_+1) + Xp(jj_,ii_+1))/4;
                yyC = (Yp(jj_,ii_) + Yp(jj_+1,ii_) + Yp(jj_+1,ii_+1) + Yp(jj_,ii_+1))/4;
                zzC = (levels(jj_,ii_,kk_)   + levels(jj_+1,ii_,kk_)   + levels(jj_+1,ii_+1,kk_)   + levels(jj_,ii_+1,kk_) + ...
                       levels(jj_,ii_,kk_+1) + levels(jj_+1,ii_,kk_+1) + levels(jj_+1,ii_+1,kk_+1) + levels(jj_,ii_+1,kk_+1))/8;
                D = sqrt((mpc(ids,1)-xxC).^2+(mpc(ids,2)-yyC).^2+(mpc(ids,3)-zzC).^2);
                [dist, ind] = mink(D, maxNUMrays);
                vox = sub2ind(sz, jj_, ii_, kk_);
                siteIDS{vox, 1} = jj_;
                siteIDS{vox, 2} = ii_;
                siteIDS{vox, 3} = kk_;
                siteIDS{vox, 4} = vox;
                siteIDS{vox, 5} = ids(ind); % Station indices affecting this voxel.
                if verbose
                    strind = join(string(ind'));
                    fprintf('voxel : x,y,z[%04i, %04i, %04i] | min dist= %8.1fm | max dist= %8.1fm | num of sites= %03i | sites id= %s\n',...
                        jj_,ii_,kk_,min(dist),max(dist),length(ind),strind); 
                end
                if ( 0 )
                    figure; hold on
                    scatter3(mpc(ids,1), mpc(ids,2), mpc(ids,3), 60, '^k', 'filled')
                    scatter3(mpc(siteIDS{vox, 5},1), mpc(siteIDS{vox, 5},2), mpc(siteIDS{vox, 5},3), 60, '^r', 'filled')
                end
            end
        end
    end
else
    if verbose, disp('centerOFVOXELSrayMINnum= 0'); end
end

if target(6)
    centerOFVOXELSrayMINnumV2 = true;
    if verbose, disp('centerOFVOXELSrayMINnumV2= 1'); end
    % Minimum rays per voxel, equivalent to the nearest stations for that voxel.
    % A voxel must be visited "maxNUMrays" times, matching "maxNUMrays" stations.
    % Used by improfile3D.m.
    Xm = nan(n,m);
    Ym = Xm;
    Zm = nan(n,m,l);
    for jj = 1 : m
        for ii = 1 : n
            Xm(ii,jj) = mean([Xp(ii,jj), Xp(ii+1,jj), Xp(ii+1,jj+1), Xp(ii,jj+1)]);
            Ym(ii,jj) = mean([Yp(ii,jj), Yp(ii+1,jj), Yp(ii+1,jj+1), Yp(ii,jj+1)]);
            for kk = 1 : l
                Zm(ii,jj,kk) = mean([levels(ii,jj,kk),   levels(ii+1,jj,kk),   levels(ii+1,jj+1,kk),   levels(ii,jj+1,kk), ...
                    levels(ii,jj,kk+1), levels(ii+1,jj,kk+1), levels(ii+1,jj+1,kk+1), levels(ii,jj+1,kk+1)]);
            end
        end
    end
    % Memory allocation.
    siteIDS = cell(n*m*l,1);
    sz = [n,m,l];
    for kk_ = 1 : l
        for ii_ = 1 : m
            for jj_ = 1 : n
                % Voxel center.
                xxC = (Xp(jj_,ii_) + Xp(jj_+1,ii_) + Xp(jj_+1,ii_+1) + Xp(jj_,ii_+1))/4;
                yyC = (Yp(jj_,ii_) + Yp(jj_+1,ii_) + Yp(jj_+1,ii_+1) + Yp(jj_,ii_+1))/4;
                zzC = (levels(jj_,ii_,kk_)   + levels(jj_+1,ii_,kk_)   + levels(jj_+1,ii_+1,kk_)   + levels(jj_,ii_+1,kk_) + ...
                       levels(jj_,ii_,kk_+1) + levels(jj_+1,ii_,kk_+1) + levels(jj_+1,ii_+1,kk_+1) + levels(jj_,ii_+1,kk_+1))/8;
                D = sqrt((mpc(ids,1)-xxC).^2+(mpc(ids,2)-yyC).^2+(mpc(ids,3)-zzC).^2);
                [dist, ind] = mink(D, maxNUMrays);
                vox = sub2ind(sz, jj_, ii_, kk_);
                siteIDS{vox, 1} = jj_;
                siteIDS{vox, 2} = ii_;
                siteIDS{vox, 3} = kk_;
                siteIDS{vox, 4} = vox;
                siteIDS{vox, 5} = ids(ind); % Station indices affecting this voxel.
                if verbose
                    strind = join(string(ind'));
                    fprintf('voxel : x,y,z[%04i, %04i, %04i] | min dist= %8.1fm | max dist= %8.1fm | num of sites= %03i | sites id= %s\n',...
                        jj_,ii_,kk_,min(dist),max(dist),length(ind),strind); 
                end
            end
        end
    end
else
    if verbose, disp('centerOFVOXELSrayMINnumV2= 0'); end
end
    
% time
tStart = tic;

num = 0;
if voxelsWALL
    num = num + n*l*2 + m*l*2 + m*n;
end
if fixAZandELcycle
    num = num + length(stepAZ)*length(stepEL);
end
if centerOFallVOXELS
    num = num + n*m*l;
end
if riseRAYSinSOMEvoxels
    L = nan(size(minLEVEL,2),1);
    for ii = 1 : size(minLEVEL,2)
        minL = maxLEVEL(1,ii) - minLEVEL(1,ii) + 1;
        L(ii) = numRAYS*(n*m*minL);
    end
    num = num + sum(L);
end

% Add zenith rays.
if addZENITH, num = num + numberOFsites; end

AZ_     = nan(num*numberOFsites,1);
EL_     = AZ_;
siteID  = AZ_; % Simplifies later calls to rays.m.
XYZsite = nan(num*numberOFsites,3);
methodID= AZ_; % Numeric method ID, used to enable or disable method weights.

k1 = 1;
k2 = 0;
for ii = 1 : numberOFsites

    siteX = mpc(ids(ii), 1);
    siteY = mpc(ids(ii), 2);
    siteZ = mpc(ids(ii), 3);

    % tmp
    AZ = nan(num,1);
    EL = AZ;
    mID= AZ;

    red = [0.6350 0.0780 0.1840];
    s = 1;
    sl = 1;
    if sum(p(1:4))>=1, figure; hold on; end 

    if voxelsWALL
        if verbose
            disp(['Preparing virtual rays (intersection with the domain voxel walls)... for site number ',num2str(ids(ii))])
        end
        % Left wall (west).
        for k = 1 : l
            for i = 1 : n
                xr = (Xp(i,1) + Xp(i+1,1))/2;
                yr = (Yp(i,1) + Yp(i+1,1))/2;
                if k == 1 && drop ~= 0
                    minl = min([levels(i,1,k+1), levels(i+1,1,k+1)]);
                    zr = minl - drop;
                else
                    zr = (levels(i,1,k) + levels(i+1,1,k) + levels(i,1,k+1) + levels(i+1,1,k+1))/4;
                end
                dh = norm([xr, yr]-[siteX, siteY]);  % Adjacent side.
                do = zr - siteZ;                     % Opposite side.
                hip= sqrt(do^2+dh^2);                % Hypotenuse.
                EL(s) = asind(do/hip);
                AZ(s) = calAz([siteX, siteY], [xr, yr]);
                mID(s)= 1;
                % Check that the computed azimuth and elevation recover the target coordinates.
                cx = xr - (siteX + dh * sind(AZ(s)));
                cy = yr - (siteY + dh * cosd(AZ(s)));
                cz = zr - (siteZ + sind(EL(s)) * hip);
                if cx > 0.001 || cy > 0.001 || cz > 0.001
                    error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                end
                if p(1)
                    if EL(s) < minEL, color = red; else, color = colors(1,:); end
                    plot3([siteX, xr], [siteY, yr], [siteZ, zr], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
        % Right wall (east).
        for k = 1 : l
            for i = 1 : n
                xr = (Xp(i,end) + Xp(i+1,end))/2;
                yr = (Yp(i,end) + Yp(i+1,end))/2;
                if k == 1 && drop ~= 0
                    minl = min([levels(i,end,k+1), levels(i+1,end,k+1)]);
                    zr = minl - drop;
                else
                    zr = (levels(i,end,k) + levels(i+1,end,k) + levels(i,end,k+1) + levels(i+1,end,k+1))/4;
                end
                dh = norm([xr, yr]-[siteX, siteY]);  % Adjacent side.
                do = zr - siteZ;                     % Opposite side.
                hip= sqrt(do^2+dh^2);                % Hypotenuse.
                EL(s) = asind(do/hip);
                AZ(s) = calAz([siteX, siteY], [xr, yr]);
                mID(s)= 1;
                % Check that the computed azimuth and elevation recover the target coordinates.
                cx = xr - (siteX + dh * sind(AZ(s)));
                cy = yr - (siteY + dh * cosd(AZ(s)));
                cz = zr - (siteZ + sind(EL(s)) * hip);
                if cx > 0.001 || cy > 0.001 || cz > 0.001
                    error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                end
                if p(1)
                    if EL(s) < minEL, color = red; else, color = colors(1,:); end
                    plot3([siteX, xr], [siteY, yr], [siteZ, zr], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
        % Upper wall (north).
        for k = 1 : l
            for i = 1 : m
                xr = (Xp(1,i) + Xp(1,i+1))/2;
                yr = (Yp(1,i) + Yp(1,i+1))/2;
                if k == 1 && drop ~= 0
                    minl = min([levels(1,i,k+1), levels(1,i+1,k+1)]);
                    zr = minl - drop;
                else
                    zr = (levels(1,i,k) + levels(1,i+1,k) + levels(1,i,k+1) + levels(1,i+1,k+1))/4;
                end
                dh = norm([xr, yr]-[siteX, siteY]);  % Adjacent side.
                do = zr - siteZ;                     % Opposite side.
                hip= sqrt(do^2+dh^2);                % Hypotenuse.
                EL(s) = asind(do/hip);
                AZ(s) = calAz([siteX, siteY], [xr, yr]);
                mID(s)= 1;
                % Check that the computed azimuth and elevation recover the target coordinates.
                cx = xr - (siteX + dh * sind(AZ(s)));
                cy = yr - (siteY + dh * cosd(AZ(s)));
                cz = zr - (siteZ + sind(EL(s)) * hip);
                if cx > 0.001 || cy > 0.001 || cz > 0.001
                    error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                end
                if p(1)
                    if EL(s) < minEL, color = red; else, color = colors(1,:); end
                    plot3([siteX, xr], [siteY, yr], [siteZ, zr], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
        % Lower wall (south).
        for k = 1 : l
            for i = 1 : m
                xr = (Xp(end,i) + Xp(end,i+1))/2;
                yr = (Yp(end,i) + Yp(end,i+1))/2;
                if k == 1 && drop ~= 0
                    minl = min([levels(end,i,k+1), levels(end,i+1,k+1)]);
                    zr = minl - drop;
                else
                    zr = (levels(end,i,k) + levels(end,i+1,k) + levels(end,i,k+1) + levels(end,i+1,k+1))/4;
                end
                dh = norm([xr, yr]-[siteX, siteY]);  % Adjacent side.
                do = zr - siteZ;                     % Opposite side.
                hip= sqrt(do^2+dh^2);                % Hypotenuse.
                EL(s) = asind(do/hip);
                AZ(s) = calAz([siteX, siteY], [xr, yr]);
                mID(s)= 1;
                % Check that the computed azimuth and elevation recover the target coordinates.
                cx = xr - (siteX + dh * sind(AZ(s)));
                cy = yr - (siteY + dh * cosd(AZ(s)));
                cz = zr - (siteZ + sind(EL(s)) * hip);
                if cx > 0.001 || cy > 0.001 || cz > 0.001
                    error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                end
                if p(1)
                    if EL(s) < minEL, color = red; else, color = colors(1,:); end
                    plot3([siteX, xr], [siteY, yr], [siteZ, zr], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
        % Domain top.
        for i = 1 : n
            for j = 1 : m
                xr = (Xp(i,j) + Xp(i+1,j) + Xp(i+1,j+1) + Xp(i,j+1))/4;
                yr = (Yp(i,j) + Yp(i+1,j) + Yp(i+1,j+1) + Yp(i,j+1))/4;
                zr = (levels(i,j,end) + levels(i+1,j,end) + levels(i+1,j+1,end) + levels(i,j+1,end))/4;
                dh = norm([xr, yr]-[siteX, siteY]);  % Adjacent side.
                do = zr - siteZ;                     % Opposite side.
                hip= sqrt(do^2+dh^2);                % Hypotenuse.
                EL(s) = asind(do/hip);
                AZ(s) = calAz([siteX, siteY], [xr, yr]);
                mID(s)= 1;
                % Check that the computed azimuth and elevation recover the target coordinates.
                cx = xr - (siteX + dh * sind(AZ(s)));
                cy = yr - (siteY + dh * cosd(AZ(s)));
                cz = zr - (siteZ + sind(EL(s)) * hip);
                if cx > 0.001 || cy > 0.001 || cz > 0.001
                    error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                end
                if p(1)
                    if EL(s) < minEL, color = red; else, color = colors(1,:); end
                    plot3([siteX, xr], [siteY, yr], [siteZ, zr], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
    end

    % Independent of voxels; use fixed AZ and EL values in a cycle.
    if fixAZandELcycle
        if verbose
            disp(['Preparing virtual rays (azimuth and elevation cycle)... for site number ',num2str(ids(ii))])
        end
        for aa = 1 : length(stepAZ)
            for ee = 1 : length(stepEL)
                AZ(s) = stepAZ(aa);
                EL(s) = stepEL(ee);
                mID(s)= 2;
                if p(2)
                    if EL(s) < minEL, color = red; else, color = colors(2,:); end
                    % Apply max domain distance divided by 2, only to inspect rays in the plot.
                    xmax = Xp(1,end)-Xp(1,1);
                    ymax = Yp(1,1)-Yp(end,1);
                    dmax = sqrt(xmax^2 + ymax^2)/2;
                    xxC = siteX + dmax*sind(AZ(s));
                    yyC = siteY + dmax*cosd(AZ(s));
                    zzC = siteZ + dmax*(1/tand(90-EL(s)));
                    plot3([siteX, xxC], [siteY, yyC], [siteZ, zzC], 'color', color, 'LineWidth', sl);
                end
                s = s + 1;
            end
        end
    end   

    % Point to the center of every voxel.
    if centerOFallVOXELS
        % In this method's plots, rays only reach voxel centers; only azimuth and elevation are needed.
        if verbose
            disp(['Preparing virtual rays (intersection with all voxels)... for site number ',num2str(ids(ii))])
        end
        for kk_ = 1 : l
            for ii_ = 1 : m
                for jj_ = 1 : n
                    % Voxel center.
                    xxC = (Xp(jj_,ii_) + Xp(jj_+1,ii_) + Xp(jj_+1,ii_+1) + Xp(jj_,ii_+1))/4;
                    yyC = (Yp(jj_,ii_) + Yp(jj_+1,ii_) + Yp(jj_+1,ii_+1) + Yp(jj_,ii_+1))/4;
                    if kk_ == 1 && drop ~= 0
                        minl = min([levels(jj_,ii_,kk_+1), levels(jj_+1,ii_,kk_+1), levels(jj_+1,ii_+1,kk_+1), levels(jj_,ii_+1,kk_+1)]);
                        zzC = minl - drop;
                    else
                        zzC = (levels(jj_,ii_,kk_)   + levels(jj_+1,ii_,kk_)   + levels(jj_+1,ii_+1,kk_)   + levels(jj_,ii_+1,kk_) + ...
                            levels(jj_,ii_,kk_+1) + levels(jj_+1,ii_,kk_+1) + levels(jj_+1,ii_+1,kk_+1) + levels(jj_,ii_+1,kk_+1))/8;
                    end
                    % Station geometry.
                    do = zzC - siteZ;                               % Opposite side.
                    dh = norm( [siteX, siteY, 0] - [xxC, yyC, 0] ); % Adjacent side.
                    hip= sqrt(do^2+dh^2);                           % Hypotenuse.
                    EL(s) = asind(do/hip);
                    AZ(s) = calAz([siteX, siteY], [xxC, yyC]);
                    mID(s)= 3;
                    % Check that the computed azimuth and elevation recover the target coordinates.
                    cx = xxC - (siteX + dh * sind(AZ(s)));
                    cy = yyC - (siteY + dh * cosd(AZ(s)));
                    cz = zzC - (siteZ + sind(EL(s)) * hip);
                    if cx > 0.001 || cy > 0.001 || cz > 0.001
                        error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                    end
                    if p(3)
                        if EL(s) < minEL, color = red; else, color = colors(3,:); end
                        plot3([siteX, xxC], [siteY, yyC], [siteZ, zzC], 'color', color, 'LineWidth', sl);
                    end
                    s = s + 1;
                end
            end
        end
    end
    
    % Increase the number of rays at preselected levels by (numRAYS-1)^2 points per voxel.
    % Target only voxels in the selected levels.
    if riseRAYSinSOMEvoxels
        if verbose
            disp(['Preparing virtual rays (rays at pre-selected levels)... for site number ',num2str(ids(ii))])
        end
        num = numRAYS; % Produces (num-1)^2 points in the voxel.
        for ll = 1 : size(minLEVEL,2)
            niv1 = minLEVEL(ll);
            niv2 = maxLEVEL(ll);
            for kk_ = niv1 : niv2
                for ii_ = 1 : m
                    for jj_ = 1 : n
                        d = norm( [Xp(jj_,ii_+1), Yp(jj_,ii_+1), 0] - [Xp(jj_,ii_), Yp(jj_,ii_), 0] )/num; 
                        v1 = levels(jj_,ii_,kk_+1)     - levels(jj_,ii_,kk_);
                        v2 = levels(jj_+1,ii_,kk_+1)   - levels(jj_+1,ii_,kk_);
                        v3 = levels(jj_+1,ii_+1,kk_+1) - levels(jj_+1,ii_+1,kk_);
                        v4 = levels(jj_,ii_+1,kk_+1)   - levels(jj_,ii_+1,kk_);
                        zC = ((v1+v2+v3+v4)/4)/num;                
                        for ss_ = 1 : num-1
                            xxC = Xp(jj_,ii_) + ss_*d;
                            yyC = Yp(jj_,ii_) - ss_*d;
                            zzC = levels(jj_,ii_,kk_) + ss_*zC; % Approximate.
                            % Station geometry.
                            do = zzC - siteZ;                               % Opposite side.
                            dh = norm( [siteX, siteY, 0] - [xxC, yyC, 0] ); % Adjacent side.
                            hip= sqrt(do^2+dh^2);                           % Hypotenuse.
                            EL(s) = asind(do/hip);
                            AZ(s) = calAz([siteX, siteY], [xxC, yyC]);
                            mID(s)= 4;
                            % Check that the computed azimuth and elevation recover the target coordinates.
                            cx = xxC - (siteX + dh * sind(AZ(s)));
                            cy = yyC - (siteY + dh * cosd(AZ(s)));
                            cz = zzC - (siteZ + sind(EL(s)) * hip);
                            if cx > 0.001 || cy > 0.001 || cz > 0.001
                                error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                            end
                            if p(4)
                                if EL(s) < minEL, color = red; else, color = colors(4,:); end
                                plot3([siteX, xxC], [siteY, yyC], [siteZ, zzC], 'color', color, 'LineWidth', sl);
                            end
                        end
                        s = s + 1;
                    end
                end
            end
        end
    end
            
    % Add zenith rays to the virtual set.
    if addZENITH
        EL(s) = 90;
        AZ(s) = 0;
        mID(s)= 0;
        siteID(s)    = ids(ii);
        XYZsite(s,1) = siteX;
        XYZsite(s,2) = siteY;
        XYZsite(s,3) = siteZ;
        s = s + 1;
    end

    % Minimum elevation relative to station height.
    idx = find(EL < minEL);
    AZ(idx) = [];
    EL(idx) = [];
    mID(idx)= [];

    % Save outputs; siteID identifies the station.
    k2 = k2 + length(AZ);
    AZ_(k1:k2) = AZ;
    EL_(k1:k2) = EL;
    methodID(k1:k2)  = mID;
    siteID(k1:k2)    = ids(ii);
    XYZsite(k1:k2,1) = siteX;
    XYZsite(k1:k2,2) = siteY;
    XYZsite(k1:k2,3) = siteZ;
    k1 = k2 + 1;
    
    if sum(p(1:4)) >= 1
        % add a 3d grid of domain
        color = 'b';
        sl = 2;
        % South-north.
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,1), levels(end,1,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,end), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        % East-west.
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,1), levels(end,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,end), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        % Verticals.
        plot3([Xp(1,1), Xp(1,1)], [Yp(1,1), Yp(1,1)], [levels(1,1,1), levels(1,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,end), Xp(1,end)], [Yp(1,end), Yp(1,end)], [levels(1,end,1), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,1)], [Yp(end,1), Yp(end,1)], [levels(end,1,1), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(end,end)], [Yp(end,end), Yp(end,end)], [levels(end,end,1), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(['SITE: ',num2str(ids(ii))])
        axis equal
    end
    
end

% Remove trailing NaN values from the first four methods.
idx = isnan(AZ_);
AZ_(idx) = [];
EL_(idx) = [];
methodID(idx) = [];
siteID(idx)   = [];
XYZsite(idx,:)= [];

% End of the first four methods.
% Methods 5 and 6 are processed outside the main loop.
if centerOFVOXELSrayMINnum
    
    num      = (n*m*l)*maxNUMrays*numberOFsites;
    AZ__     = nan(num,1);
    EL__     = AZ__;
    siteID_  = AZ__; % Simplifies later calls to rays.m.
    XYZsite_ = nan(num,3);
    methodID_= AZ__; % Numeric method ID, used to enable or disable method weights.
    
    k3 = 1;
    sl = 1;
    if p(5), figure; hold on; title('centerOFVOXELSrayMINnum'); end
    
    % Nearest stations for each voxel were determined at the beginning.
    for kk_ = 1 : l
        for ii_ = 1 : m
            for jj_ = 1 : n
                % Target voxel center.
                xxC = (Xp(jj_,ii_) + Xp(jj_+1,ii_) + Xp(jj_+1,ii_+1) + Xp(jj_,ii_+1))/4;
                yyC = (Yp(jj_,ii_) + Yp(jj_+1,ii_) + Yp(jj_+1,ii_+1) + Yp(jj_,ii_+1))/4;
                if kk_ == 1 && drop ~= 0
                    minl = min([levels(jj_,ii_,kk_+1), levels(jj_+1,ii_,kk_+1), levels(jj_+1,ii_+1,kk_+1), levels(jj_,ii_+1,kk_+1)]);
                    zzC = minl - drop;
                else
                    zzC = (levels(jj_,ii_,kk_)   + levels(jj_+1,ii_,kk_)   + levels(jj_+1,ii_+1,kk_)   + levels(jj_,ii_+1,kk_) + ...
                        levels(jj_,ii_,kk_+1) + levels(jj_+1,ii_,kk_+1) + levels(jj_+1,ii_+1,kk_+1) + levels(jj_,ii_+1,kk_+1))/8;
                end
                
                vox = sub2ind([n,m,l], jj_, ii_, kk_);
                for ii = 1 : maxNUMrays
                    % Station loop.
                    site  = siteIDS{vox, 5}(ii);
                    siteX = mpc(siteIDS{vox, 5}(ii), 1);
                    siteY = mpc(siteIDS{vox, 5}(ii), 2);
                    siteZ = mpc(siteIDS{vox, 5}(ii), 3);
                    % Station "ii" influences voxel "vox".
                    % Add the ray.
                    do = zzC - siteZ;                               % Opposite side.
                    dh = norm( [siteX, siteY, 0] - [xxC, yyC, 0] ); % Adjacent side.
                    hip= sqrt(do^2+dh^2);                           % Hypotenuse.
                    EL = asind(do/hip);
                    AZ = calAz([siteX, siteY], [xxC, yyC]);
                    % Check that the computed azimuth and elevation recover the target coordinates.
                    cx = xxC - (siteX + dh * sind(AZ));
                    cy = yyC - (siteY + dh * cosd(AZ));
                    cz = zzC - (siteZ + sind(EL) * hip);
                    if cx > 0.001 || cy > 0.001 || cz > 0.001
                        error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                    end
                    if p(5)
                        if EL < minEL, color = red; else, color = colors(5,:); end
                        plot3([siteX, xxC], [siteY, yyC], [siteZ, zzC], 'color', color, 'LineWidth', sl);
                    end
                    % Save outputs; siteID identifies the station.
                    if EL >= minEL
                        AZ__(k3) = AZ;
                        EL__(k3) = EL;
                        methodID_(k3)  = 5;
                        siteID_(k3)    = site;
                        XYZsite_(k3,1) = siteX;
                        XYZsite_(k3,2) = siteY;
                        XYZsite_(k3,3) = siteZ;
                        k3 = k3 + 1;
                    end
                end
            end
        end
    end
    
    % Add virtual zenith rays; one zenith per station is enough.
    if addZENITH
        for ii = 1 : numberOFsites
            EL__(k3)      = 90;
            AZ__(k3)      = 0;
            methodID_(k3) = 0;
            siteID_(k3)   = ids(ii);
            XYZsite_(k3,1)= mpc(ids(ii),1);
            XYZsite_(k3,2)= mpc(ids(ii),2);
            XYZsite_(k3,3)= mpc(ids(ii),3);
            k3 = k3 + 1;
        end
    end
    
    % Remove trailing NaN values for this method.
    idx = isnan(AZ__);
    AZ__(idx)      = [];
    EL__(idx)      = [];
    methodID_(idx) = [];
    siteID_(idx)   = [];
    XYZsite_(idx,:)= [];
    
    % Merge methods.
    AZ_     = cat(1, AZ_, AZ__);
    EL_     = cat(1, EL_, EL__);
    siteID  = cat(1, siteID, siteID_);
    XYZsite = cat(1, XYZsite, XYZsite_);
    methodID= cat(1, methodID, methodID_);
    
    if p(5)
        % add a 3d grid of domain
        color = 'b';
        sl = 2;
        % South-north.
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,1), levels(end,1,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,end), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        % East-west.
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,1), levels(end,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,end), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        % Verticals.
        plot3([Xp(1,1), Xp(1,1)], [Yp(1,1), Yp(1,1)], [levels(1,1,1), levels(1,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,end), Xp(1,end)], [Yp(1,end), Yp(1,end)], [levels(1,end,1), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,1)], [Yp(end,1), Yp(end,1)], [levels(end,1,1), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(end,end)], [Yp(end,end), Yp(end,end)], [levels(end,end,1), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        % Stations.
        scatter3(mpc(ids,1), mpc(ids,2), mpc(ids,3), 150, '^k', 'filled')
        axis equal
    end
    
end

% Methods 5 and 6 are processed outside the main loop.
if centerOFVOXELSrayMINnumV2
    
    num      = (n*m*l)*maxNUMrays*numberOFsites;
    AZ__     = nan(num,1);
    EL__     = AZ__;
    siteID_  = AZ__; % Simplifies later calls to rays.m.
    XYZsite_ = nan(num,3);
    methodID_= AZ__; % Numeric method ID, used to enable or disable method weights.
   
    if ~centerOFVOXELSrayMINnum
        k3 = 1;
    end
    sl = 1;
    if p(6), figure; hold on; title('centerOFVOXELSrayMINnumV2'); end
    
    voxNUMrays = zeros(n*m*l,1);
    %voxNUMrays_ELmin = voxNUMrays; % For testing.
    for kk_ = 1 : l
        for ii_ = 1 : m
            for jj_ = 1 : n
                % Target voxel center.
                xxC = (Xp(jj_,ii_) + Xp(jj_+1,ii_) + Xp(jj_+1,ii_+1) + Xp(jj_,ii_+1))/4;
                yyC = (Yp(jj_,ii_) + Yp(jj_+1,ii_) + Yp(jj_+1,ii_+1) + Yp(jj_,ii_+1))/4;
                if kk_ == 1 && drop ~= 0
                    minl = min([levels(jj_,ii_,kk_+1), levels(jj_+1,ii_,kk_+1), levels(jj_+1,ii_+1,kk_+1), levels(jj_,ii_+1,kk_+1)]);
                    zzC = minl - drop;
                else
                    zzC = (levels(jj_,ii_,kk_)   + levels(jj_+1,ii_,kk_)   + levels(jj_+1,ii_+1,kk_)   + levels(jj_,ii_+1,kk_) + ...
                        levels(jj_,ii_,kk_+1) + levels(jj_+1,ii_,kk_+1) + levels(jj_+1,ii_+1,kk_+1) + levels(jj_,ii_+1,kk_+1))/8;
                end                
                
                vox = sub2ind([n,m,l], jj_, ii_, kk_);
                for ii = 1 : maxNUMrays
                    if voxNUMrays(vox) < maxNUMrays
                        % Station loop.
                        site  = siteIDS{vox, 5}(ii);
                        siteX = mpc(siteIDS{vox, 5}(ii), 1);
                        siteY = mpc(siteIDS{vox, 5}(ii), 2);
                        siteZ = mpc(siteIDS{vox, 5}(ii), 3);
                        % Station "ii" influences voxel "vox".
                        % Add the ray.
                        do = zzC - siteZ;                               % Opposite side.
                        dh = norm( [siteX, siteY, 0] - [xxC, yyC, 0] ); % Adjacent side.
                        hip= sqrt(do^2+dh^2);                           % Hypotenuse.
                        EL = asind(do/hip);
                        AZ = calAz([siteX, siteY], [xxC, yyC]);
                        % Check that the computed azimuth and elevation recover the target coordinates.
                        cx = xxC - (siteX + dh * sind(AZ));
                        cy = yyC - (siteY + dh * cosd(AZ));
                        cz = zzC - (siteZ + sind(EL) * hip);
                        if cx > 0.001 || cy > 0.001 || cz > 0.001
                            error(['Check the azimuth and elevation for site ',num2str(ids(ii))])
                        end
                        if p(6)
                            if EL < minEL, color = red; else, color = colors(6,:); end
                            plot3([siteX, xxC], [siteY, yyC], [siteZ, zzC], 'color', color, 'LineWidth', sl);
                        end                     
                        % Extract all voxels crossed by the ray.
                        [~, isub] = improfile3D(Xm, Ym, Zm, [siteX, siteY, siteZ], [xxC, yyC, zzC]);
                        voxNUMrays(isub) = voxNUMrays(isub) + 1;
                        %
                        % Save outputs; siteID identifies the station.
                        if EL >= minEL
                            %voxNUMrays_ELmin(vox) = voxNUMrays_ELmin(vox) + 1;
                            AZ__(k3) = AZ;
                            EL__(k3) = EL;   
                            methodID_(k3)  = 6;
                            siteID_(k3)    = site;
                            XYZsite_(k3,1) = siteX;
                            XYZsite_(k3,2) = siteY;
                            XYZsite_(k3,3) = siteZ;
                            k3 = k3 + 1;
                        end
                    end
                end
            end
        end
    end
    
    % Add virtual zenith rays; one zenith per station is enough.
    if addZENITH
        for ii = 1 : numberOFsites
            EL__(k3)      = 90;
            AZ__(k3)      = 0;
            methodID_(k3) = 0;
            siteID_(k3)   = ids(ii);
            XYZsite_(k3,1)= mpc(ids(ii),1);
            XYZsite_(k3,2)= mpc(ids(ii),2);
            XYZsite_(k3,3)= mpc(ids(ii),3);
            k3 = k3 + 1;
        end
    end
    
    % Remove trailing NaN values for this method.
    idx = isnan(AZ__);
    AZ__(idx)      = [];
    EL__(idx)      = [];
    methodID_(idx) = [];
    siteID_(idx)   = [];
    XYZsite_(idx,:)= [];
    
    % Merge methods.
    AZ_     = cat(1, AZ_, AZ__);
    EL_     = cat(1, EL_, EL__);
    siteID  = cat(1, siteID, siteID_);
    XYZsite = cat(1, XYZsite, XYZsite_);
    methodID= cat(1, methodID, methodID_);
        
    if p(6)
        % add a 3d grid of domain
        color = 'b';
        sl = 2;
        % South-north.
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,1), levels(end,1,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(end,1)], [Yp(1,1), Yp(end,1)], [levels(1,1,end), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(1,end)], [Yp(end,end), Yp(1,end)], [levels(end,end,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        % East-west.
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,1), levels(1,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,1), Xp(1,end)], [Yp(1,1), Yp(1,end)], [levels(1,1,end), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,1), levels(end,end,1)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,end)], [Yp(end,1), Yp(end,end)], [levels(end,1,end), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        % Verticals.
        plot3([Xp(1,1), Xp(1,1)], [Yp(1,1), Yp(1,1)], [levels(1,1,1), levels(1,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(1,end), Xp(1,end)], [Yp(1,end), Yp(1,end)], [levels(1,end,1), levels(1,end,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,1), Xp(end,1)], [Yp(end,1), Yp(end,1)], [levels(end,1,1), levels(end,1,end)], 'color', color, 'LineWidth', sl);
        plot3([Xp(end,end), Xp(end,end)], [Yp(end,end), Yp(end,end)], [levels(end,end,1), levels(end,end,end)], 'color', color, 'LineWidth', sl);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        % Stations.
        scatter3(mpc(ids,1), mpc(ids,2), mpc(ids,3), 150, '^k', 'filled')
        axis equal
    end
    
end

% Sort only to improve visualization of matrix W when it is used.
[AZ_, ij] = sort(AZ_);
EL_= EL_(ij);
siteID = siteID(ij);
XYZsite = XYZsite(ij,:);
methodID = methodID(ij);

% Terminal output.
tEnd = toc(tStart);
if verbose
    disp(['... done in ',num2str(tEnd, '%.3f'),' sec '])
end

end


function azv = calAz(ipoint, fpoint)
% calAz computes azimuth from cartographic north between two points.

% x
x1 = ipoint(1);
y1 = ipoint(2);
% y
x2 = fpoint(1);
y2 = fpoint(2);

u = (x2 - x1);
v = (y2 - y1);
if u >= 0 && v >= 0, c = 0; end   % 1Q
if u >  0 && v <  0, c = 180; end % 2Q
if u <  0 && v <  0, c = 180; end % 3Q
if u <  0 && v >  0, c = 360; end % 4Q
azv = atand(u/v) + c;
end

function [C, point_list] = improfile3D(Xm,Ym,Zm,p1,p2)
% for testing
% p1=[siteX, siteY, siteZ]; p2=[xxC, yyC, zzC];

n = size(Xm, 1); % Rows.
m = size(Xm, 2); % Columns.
l = size(Zm, 3);

[xx,yy] = meshgrid(1:m,1:n);
zz = 1:l;

%figure; hold on
%scatter(Xm(:), Ym(:), 50)
%scatter(p1(1), p1(2), 50, 'red')

F = scatteredInterpolant(Xm(:), Ym(:), xx(:), 'nearest', 'nearest');
point1x = F(p1(1), p1(2));
point2x = F(p2(1), p2(2));
F = scatteredInterpolant(Xm(:), Ym(:), yy(:), 'nearest', 'nearest');
point1y = F(p1(1), p1(2));
point2y = F(p2(1), p2(2));
v = squeeze(Zm(point1y, point1x, :));
point1z = interp1(v, zz, p1(3), 'nearest', 'extrap');
v = squeeze(Zm(point2y, point2x, :));
point2z = interp1(v, zz, p2(3), 'nearest', 'extrap');
p1 = [point1y, point1x, point1z];
p2 = [point2y, point2x, point2z];

% Euclidian distance
dist_euc = norm(p1 - p2);
if dist_euc == 0
    % The ray points to the center of the voxel containing the station.
    point_list = sub2ind([n, m, l], p1(1), p1(2), p1(3));
    C = p1;
else
    % Number of intervals.
    n_intervalle = round(dist_euc*12);
    step = (p1 - p2)/n_intervalle;
    pix_coords = zeros(n_intervalle, 3);
    for cp = 0:n_intervalle
        pix_coords(cp+1, :) = round(p2 + cp*step);
    end
    point_list = sub2ind([n, m, l],pix_coords(:,1),pix_coords(:,2),pix_coords(:,3));
    point_list = unique(point_list);
    
    % Point coordinates
    x = pix_coords(:,1);
    y = pix_coords(:,2);
    z = pix_coords(:,3);
    
    C = unique([x, y, z],'rows');
end
end

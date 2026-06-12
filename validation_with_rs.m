function [rmse_, bias_, corr_, rhov_tomo, rs_rhov, hgt_tomo] = validation_with_rs(Xp, Yp, levels, dt, tomoX, path2rs, rs_coord, p)
% validation_with_rs compares tomographic water-vapor density with radiosonde data.
% Radiosonde profiles are voxelized by rs_voxalization_1d.m before computing
% RMSE, bias, and correlation along the sonde path.
% Positive bias means TOMO has higher rho_v than the radiosonde.
%
% rhov_tomo, rs_rhov, hgt_tomo : outputs used for plots in the main program.
%
% for testing
% dt=rdtimes(t);tomoX=TOMOdata.x{t};tomoCX=TOMOdata.cx{t};path2rs=FIXdata.paths{3};rs_coord=[xOBS(t),yOBS(t)];
% 

rmse_ = nan; bias_ = nan; corr_ = nan; rhov_tomo = nan; hgt_tomo = nan;

% Coordinates in meters.
rsX = rs_coord(1); % Sonde coordinate.
rsY = rs_coord(2);

% Voxelize the radiosonde.
[rs_rhov, ~, rs_voxels, ~] = rs_voxalization_1d(path2rs, dt, [rsX,rsY], Xp, Yp, levels, 0);

if sum(rs_rhov, [], 'omitnan') > 0
    
    n = size(Xp,1)-1;
    m = size(Xp,2)-1;
    l = size(levels,3)-1;
    
    % Reconstruct the 3D matrix.
    tomoSOLUTION = nan(n,m,l);    
    v = 1;
    for z_ = 1 : l
        for x_ = 1 : m
            for y_ = 1 : n
                tomoSOLUTION(y_,x_,z_) = tomoX(v);              
                v = v + 1;
            end
        end
    end
    
    % Compare with TOMO at the sonde location.
    rhov_tomo = nan(size(rs_voxels,1),1);
    hgt_tomo  = rhov_tomo;
    for j = 1 : size(rs_voxels,1)
        rhov_tomo(j) = tomoSOLUTION(rs_voxels(j,1),rs_voxels(j,2),rs_voxels(j,3));        
        hgt_tomo(j) = levels(rs_voxels(j,1),rs_voxels(j,2),rs_voxels(j,3));
    end
    
    id0 = find(isnan(rs_rhov));
    if ~isempty(id0)
        id1 = 1:id0(1)-1;
    else
        id1 = 1:length(rs_rhov);
    end
    rmse_ = sqrt(mean((rhov_tomo(id1) - rs_rhov(id1)).^2))*1000; % g/m^3
    bias_ = (mean(rhov_tomo(id1)) - mean(rs_rhov(id1)))*1000;
    corr_ = corr(rhov_tomo(id1), rs_rhov(id1));
    
    if p
        figure('color', 'w'); hold on
        plot(rhov_tomo.*1000, hgt_tomo);
        plot(rs_rhov.*1000, hgt_tomo);
        legend('TOMO', 'RS')
        xlabel('rho_v (g/m^3)'); ylabel('Height (m)')
    end
    
end
return

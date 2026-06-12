% plot_hovmoller_tomo_rs
% Build a two-panel Hovmoller diagram from TOMO and radiosonde 1D profiles.
% The first panel shows TOMOdata.tomo_rs_rhov_1d and the second panel shows
% TOMOdata.rs_rhov_1d, using TOMOdata.rs_hgt_1d as the vertical coordinate.

clear; close all; clc

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

matFile = fullfile(scriptDir, 'solutions', 'hong_kong', ...
    'hong_kong_minEL03_a60min_x09y06z10_rx6km_ry6km_hvRS_rot25_v1.mat');

rhoScale = 1000;           % kg/m^3 to g/m^3.
rhoUnits = 'g m^{-3}';
heightScale = 1000;        % m to km.
heightUnits = 'km';
saveFigure = false;
outFigure = fullfile(scriptDir, 'solutions', 'hong_kong', ...
    'hovmoller_tomo_rs_hong_kong_minEL03_a60min_x09y06z10_rx6km_ry6km_hvRS_rot25_v1.png');

S = load(matFile, 'TOMOdata', 'FIXdata');

requiredFields = {'tomo_rs_rhov_1d', 'rs_rhov_1d', 'rs_hgt_1d'};
for k = 1:numel(requiredFields)
    if ~isfield(S.TOMOdata, requiredFields{k})
        error('Missing TOMOdata.%s in %s', requiredFields{k}, matFile)
    end
end

tomoRhov = profilesToMatrix(S.TOMOdata.tomo_rs_rhov_1d) .* rhoScale;
rsRhov = profilesToMatrix(S.TOMOdata.rs_rhov_1d) .* rhoScale;
rsHgt = profilesToMatrix(S.TOMOdata.rs_hgt_1d);

nTimes = size(tomoRhov, 2);
if size(rsRhov, 2) ~= nTimes || size(rsHgt, 2) ~= nTimes
    error('TOMO, RS, and height profiles must have the same number of time steps.')
end

heightM = representativeHeight(rsHgt);
height = heightM ./ heightScale;
timeNum = profileTimes(S, nTimes);

validRows = ~(all(isnan(tomoRhov), 2) & all(isnan(rsRhov), 2) & isnan(heightM));
tomoRhov = tomoRhov(validRows, :);
rsRhov = rsRhov(validRows, :);
height = height(validRows);

clim = commonColorLimits(tomoRhov, rsRhov);

figure('Color', 'w', 'Position', [100, 100, 1200, 720])
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1 = nexttile;
plotHovmollerPanel(timeNum, height, tomoRhov, ...
    'Tomography at radiosonde column', heightUnits, rhoUnits, clim)

ax2 = nexttile;
plotHovmollerPanel(timeNum, height, rsRhov, ...
    'Radiosonde profile', heightUnits, rhoUnits, clim)
xlabel(ax2, 'Time')

linkaxes([ax1, ax2], 'x')

if saveFigure
    exportgraphics(gcf, outFigure, 'Resolution', 200)
    disp(['Saved figure: ', outFigure])
end

function M = profilesToMatrix(profiles)
% Convert a cell array of vertical profiles into a levels-by-time matrix.
if ~iscell(profiles)
    M = double(profiles);
    if isvector(M)
        M = M(:);
    end
    return
end

nTimes = numel(profiles);
nLevels = max(cellfun(@numel, profiles));
M = nan(nLevels, nTimes);
for t = 1:nTimes
    if isempty(profiles{t})
        continue
    end
    v = double(profiles{t}(:));
    M(1:numel(v), t) = v;
end
end

function heightM = representativeHeight(rsHgt)
% Use a single vertical coordinate; if heights vary with time, average them.
validCols = find(any(~isnan(rsHgt), 1));
if isempty(validCols)
    error('TOMOdata.rs_hgt_1d does not contain valid height values.')
end

ref = rsHgt(:, validCols(1));
sameGrid = true;
for k = validCols
    idx = ~isnan(ref) & ~isnan(rsHgt(:, k));
    if any(abs(rsHgt(idx, k) - ref(idx)) > 1e-6)
        sameGrid = false;
        break
    end
end

if sameGrid
    heightM = ref;
else
    heightM = mean(rsHgt, 2, 'omitnan');
    disp('Vertical levels vary with time; using the time-mean height for each level.')
end
end

function timeNum = profileTimes(S, nTimes)
% Prefer the center of FIXdata.times intervals; otherwise use profile index.
if isfield(S, 'FIXdata') && isfield(S.FIXdata, 'times') && size(S.FIXdata.times, 1) == nTimes
    timeNum = S.FIXdata.times(:, 2);
elseif isfield(S, 'FIXdata') && isfield(S.FIXdata, 'rtimes') && numel(S.FIXdata.rtimes) == nTimes
    timeNum = S.FIXdata.rtimes(:);
else
    timeNum = (1:nTimes)';
end
end

function clim = commonColorLimits(A, B)
% Use shared limits so both panels are directly comparable.
v = [A(:); B(:)];
v = v(isfinite(v));
if isempty(v)
    clim = [0, 1];
elseif min(v) == max(v)
    clim = min(v) + [-0.5, 0.5];
else
    clim = [min(v), max(v)];
end
end

function plotHovmollerPanel(timeNum, height, data, panelTitle, heightUnits, rhoUnits, clim)
contourf(timeNum, height, data, 24, 'LineStyle', 'none')
set(gca, 'YDir', 'normal')
ylim([min(height, [], 'omitnan'), max(height, [], 'omitnan')])
caxis(clim)
colormap(turbo)
grid on
box on
title(panelTitle)
ylabel(['Height (', heightUnits, ')'])
cb = colorbar;
ylabel(cb, ['\rho_v (', rhoUnits, ')'])

if numel(timeNum) > 1 && all(timeNum > 700000)
    datetick('x', 'dd-mmm HH:MM', 'keeplimits')
end
end

function [x, y] = gauss(lon0, lat0, lon, lat)
% gauss converts geographic coordinates to local Gauss-Kruger transverse Mercator coordinates.
% Inputs and outputs are vectorized; returned coordinates are in meters.

% ------------------------------------------------------------
% WGS84
% ------------------------------------------------------------
a = 6378.137;
f = 0.0033528106647474805;

equad = f*(2 - f);
b = a*sqrt(1 - equad);

% ------------------------------------------------------------
% Unit conversions.
% ------------------------------------------------------------
lat0 = deg2rad(lat0);
lon0 = deg2rad(lon0);

lat = deg2rad(lat);
lon = deg2rad(lon);

lon = lon - lon0;

% ------------------------------------------------------------
% Auxiliary terms.
% ------------------------------------------------------------
sin_lat = sin(lat);
cos_lat = cos(lat);
tan_lat = tan(lat);

N  = a ./ sqrt(1 - equad .* (sin_lat.^2));
RO = a*(1 - equad) ./ (1 - equad*(sin_lat.^2)).^(3/2);

% Coefficients.
k1 = (N./RO) + 4*(N.^2)./(RO.^2) - (tan_lat.^2);
k2 = (N./RO) - (tan_lat.^2);
k3 = (N./RO).*(14 - 58*(tan_lat.^2)) + ...
     40*(tan_lat.^2) + (tan_lat.^4) - 9;

% ------------------------------------------------------------
% Meridian arc.
% ------------------------------------------------------------
n = (a - b) / (a + b);

a0 = 1 + (n^2)/4 + (n^4)/64;
a2 = (3/2)*(n - n^3/8);
a4 = (15/16)*(n^2 - n^4/4);
a6 = (35/48)*n^3;

s1 = a/(1+n)*(a0*lat0 - a2*sin(2*lat0) + a4*sin(4*lat0) - a6*sin(6*lat0));

s2 = a/(1+n)*(a0*lat ...
     - a2*sin(2*lat) ...
     + a4*sin(4*lat) ...
     - a6*sin(6*lat));

s = s2 - s1;

% ------------------------------------------------------------
% Projection.
% ------------------------------------------------------------
x = lon .* N .* cos_lat ...
  + (lon.^3)/6 .* N .* (cos_lat.^3) .* k2 ...
  + (lon.^5)/120 .* N .* (cos_lat.^5) .* k3;

y = s ...
  + (lon.^2)/2 .* N .* sin_lat .* cos_lat ...
  + (lon.^4)/24 .* N .* sin_lat .* (cos_lat.^3) .* k1;

% ------------------------------------------------------------
% Convert kilometers to meters.
% ------------------------------------------------------------
x = x * 1000;
y = y * 1000;

end

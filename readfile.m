%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% example to read a netCDF file
% mer_his_030.nc is one file for one year
%   with monthly output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% useful path
% path to the model output
fname = 'C:\Users\61414\Downloads\mer_his_0030.nc';
% path to the model grid file
gname = 'C:\Users\61414\Downloads\misom020_grd.nc';

% load bathymetry for the continental shelf around the Mertz
% check the help on ncread to understand better how it works
% but the numbers between [...] are the indexes of a region, 
% you can check with panoply what these indexes are, as we cannot just 
% enter the lat/lon
bathy = ncread(gname,'h',[40 1],[180 140]); 
lat = ncread(gname,'lat_rho',[40 1],[180 140]); 
lon = ncread(gname,'lon_rho',[40 1],[180 140]); 
% load ice draft for the same area
ice_draft = ncread(gname,'zice',[40 1],[180 140]);
% loading masks
land_mask = ncread(gname,'mask_rho',[40 1],[180 140]);
ice_mask = ncread(gname,'mask_zice',[40 1],[180 140]);
% load surface potential temperature
theta_surf = ncread(fname,'temp',[40 1 31 1],[180 140 1 Inf]);
% load temperature at 1 location for the full water column
theta_loc = ncread(fname,'temp',[156 83 1 1],[1 1 Inf Inf]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how to find the depth
% with the terrain following 
% I didn't have time to write something here!
% But next time, play with the above, plot some figures
% just to practice

hhhhhhh

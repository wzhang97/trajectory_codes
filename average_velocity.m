%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Grounded model with average velocity
%
%    Author: weiyan ZHANG
%    Created: Oct 2018
%    Updated:
%    Last Modif by:
%    Last Modif on: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
addpath(genpath('/v_mertz/iceberg_project/ROMS_data/'));

% load grid information from the model
gname = '/v_mertz/iceberg_project/ROMS_data/misom020_grd_small.nc';
uname = '/v_mertz/iceberg_project/ROMS_data/uvel_small.nc';
vname = '/v_mertz/iceberg_project/ROMS_data/vvel_small.nc';
lat_rho = ncread(gname,'lat_rho'); % 96 * 120
lon_rho = ncread(gname,'lon_rho');
x_rho = ncread(gname,'x_rho');
y_rho = ncread(gname,'y_rho');
[dx_rho,dy_rho,dz_rho] = calc_rhogrid_dxdydz(gname); % calcualte the length of each edges of a cell
depths_rho = calc_depths_rho(gname); % depths of each level (96 * 120 * 31)
mask_zice = ncread(gname,'mask_zice'); % mask on ice-shelf rho-points (option(0) = 'land', option(1) = 'water')
mask_land = ncread(gname,'mask_rho'); % mask on RHO-points (option(0) = 'land', option(1) = 'water')
uo_read = ncread(uname,'u'); % meter second-1 (95x119x31x1456)
vo_read = ncread(vname,'v');
sustr_read = ncread(uname,'sustr'); % newton meter-2 (95x119x1456) 
svstr_read = ncread(vname,'svstr');
[X,Y] = meshgrid([1:size(lat_rho,2)],[1:size(lon_rho,1)]);
X_lim = x_rho(96,120);
Y_lim = y_rho(96,120);

% indicate the initial location of the particle
x_ini = 60;
y_ini = 80;
x = x_ini;
y = y_ini;
xx = x_rho(x,y);
yy = y_rho(x,y);

% set variable useful for the trajectory equation
U = 0;
V = 0;

% ocean coeff
Co = 0.85; % dimensionless coefficient of resistance
Cdo_skin = 0.0055; % ocean-iceberg friction coeff
rho_h2o = 1028; % seawater density kg m-3

% atmospheric coeff
rho_air = 1.225; % air density kg m-3
Ca = 0.4; % dimensionless coefficient of resistance
Cda_skin = 0.0022; % air-iceberg friction coeff
Cd = 1.25e-3;  % dimensionless air-sea friction coeff

% coriolis coeff
w_r = 7.2921e-5;
f = 2 * w_r * sin(-67 * pi / 180);
    
% size of the iceberg 
rho_icb = 916.7; % kg m-3
depth = 10; % m
depth_icb_under = depth * (rho_icb / rho_h2o);
depth_icb_above = depth - depth_icb_under;
length = 5000; % assume width is equal to length
dA_o = length * depth_icb_under; % m^2 ideally per face 
dA_a = length * depth_icb_above; % m^2 ideally per face
Ad = length ^ 2;
M = rho_icb * length ^ 2 * depth; % kg

% data set
step = 500;
time_all = zeros(1,step);
U_all = zeros(1,step);
V_all = zeros(1,step);
x_all = zeros(1,step);
y_all = zeros(1,step);
xx_all = zeros(1,step);
yy_all = zeros(1,step);
Fa_all = zeros(1,step);
Fo_all = zeros(1,step);
Fcp_all = zeros(1,step);
vel_all = zeros(1,step);
z_all = zeros(1,step);
uo_cst_all = zeros(1,step);
vo_cst_all = zeros(1,step);
ua_all = zeros(1,step);
va_all = zeros(1,step);
uo_cst_surf_all = zeros(1,step);
vo_cst_surf_all = zeros(1,step);
s_y = zeros(96,120);
s_x = zeros(96,120);
area_cell_ocean_u = zeros(96,120,31);
area_cell_ocean_v = zeros(96,120,31);
area_cell_bottom = zeros(96,120);
area_cell_air_u = zeros(96,120);
area_cell_air_v = zeros(96,120);

% timestep
dt = 6 * 60 * 60; % s

% create a new matrix of ocean velocity(94*119)
for XI = 1:94
    for ETA = 1:119
        uo_cst_rho(XI,ETA,:,:) = (uo_read(XI,ETA,:,:) + uo_read(XI+1,ETA,:,:)) / 2;
        vo_cst_rho(XI,ETA,:,:) = (vo_read(XI,ETA,:,:) + vo_read(XI+1,ETA,:,:)) / 2;
    end
end

% do a loop
for i = 1:step

% position of iceberg
    tmp_z = depths_rho(x,y,:);
    [closest_diff_z, index_z] = min(abs((-tmp_z) - depth_icb_under)); % closest_diff is the difference between the bottom and the closest layer, index is the number of the closest layer
    icb_x = 0;
    icb_y = 0;
    for index_x = x:95        
        icb_x = dx_rho(index_x,y) + icb_x;
        if length < icb_x
           break
        end
    end
    for index_y = y:120
        icb_y = dy_rho(index_x,index_y) + icb_y;
        if length < icb_y
           break
        end
    end
    s_x(index_x,index_y) = icb_x;
    s_y(index_x,index_y) = icb_y;
% area percentage under the ocean at u component(side)
    for b = y:index_y - 1
        for c = 1:index_z
            area_cell_ocean_u(x,b,c) = dy_rho(x,b) * dz_rho(x,b,c);
        end
    end
    area_total_ocean_u_1 = sum(area_cell_ocean_u(x,:,:));
    area_total_ocean_u = sum(area_total_ocean_u_1);

% area percentage under the ocean at v component(side)
    for a2 = x:index_x - 1
        for c2 = 1:index_z
            area_cell_ocean_v(a2,y,c2) = dy_rho(a2,y) * dz_rho(a2,y,c2);
        end
    end
    area_total_ocean_v_1 = sum(area_cell_ocean_v(:,y,:));
    area_total_ocean_v = sum(area_total_ocean_v_1);

% average ocean velocity on one side
    uo_average = 0; % meter second-1
    vo_average = 0;
    for a3 = x:index_x - 1
        for b3 = y:index_y - 1
            for c3 = 1:index_z                
                uo_ini = uo_cst_rho(a3,b3,c3,i) * (area_cell_ocean_u(a3,b3,c3) / area_total_ocean_u);
                vo_ini = vo_cst_rho(a3,b3,c3,i) * (area_cell_ocean_v(a3,b3,c3) / area_total_ocean_v);
                uo_average = uo_average + uo_ini;
                vo_average = vo_average + vo_ini;
            end
        end
    end

    uo_cst_all(i) = uo_average;
    vo_cst_all(i) = vo_average;

% area percentage at the bottom
    for a4 = x:index_x - 1
        for b4 = y:index_y - 1
            area_cell_bottom(a4,b4) = dx_rho(a4,b4) * dy_rho(a4,b4);

        end
    end
    area_total_bottom_1 = sum(area_cell_bottom);
    area_total_bottom = sum(area_total_bottom_1);

% average ocean velocity at bottom
    uo_bottom_average = 0; % meter second-1
    vo_bottom_average = 0;
    for a5 = x:index_x - 1
        for b5 = y:index_y - 1
            uo_bottom_ini = uo_cst_rho(a5,b5,index_z,i) * (area_cell_bottom(a5,b5) / area_total_bottom);
            vo_bottom_ini = vo_cst_rho(a5,b5,index_z,i) * (area_cell_bottom(a5,b5) / area_total_bottom);
            uo_bottom_average = uo_bottom_average + uo_bottom_ini;
            vo_bottom_average = vo_bottom_average + vo_bottom_ini;
        end
    end

    uo_skin_all(i) = uo_bottom_average;
    vo_skin_all(i) = vo_bottom_average;

% ua = sqrt(sustr(x,y,step) / (rho_air * Cd)); m s-2

% area percentage above the ocean at u component
    for b6 = y:index_y - 1
        area_cell_air_u(x,b6) = dy_rho(x,b6) * depth_icb_above;
    end
    area_total_air_u = sum(area_cell_air_u(x,:));

% area percentage above the ocean at v component
    for a7 = x:index_x - 1
        area_cell_air_v(a7,y) = dy_rho(a7,y) * depth_icb_above;
    end
    area_total_air_v = sum(area_cell_air_v(:,y));

% average atmospheric velocity(side)
    ua_average = 0; % meter second-1
    va_average = 0;
    for a8 = x:index_x - 1
        for b8 = y:index_y - 1
            ua_ini = sqrt(abs(sustr_read(a8,b8,i)) / (rho_air * Cd)) * (area_cell_air_u(a8,b8) / area_total_air_u);
            va_ini = sqrt(abs(sustr_read(a8,b8,i)) / (rho_air * Cd)) * (area_cell_air_v(a8,b8) / area_total_air_v);
            ua_average = ua_average + ua_ini;
            va_average = va_average + va_ini;
        end
    end
    ua_all(i) = ua_average;
    va_all(i) = va_average;

% average ocean velocity at top
    ua_top_average = 0; % meter second-1
    va_top_average = 0;
    for a9 = x:index_x - 1
        for b9 = y:index_y - 1
            ua_top_ini = sqrt(abs(sustr_read(a9,b9,i)) / (rho_air * Cd)) * (area_cell_bottom(a9,b9) / area_total_bottom);
            va_top_ini = sqrt(abs(sustr_read(a9,b9,i)) / (rho_air * Cd)) * (area_cell_bottom(a9,b9) / area_total_bottom);
            ua_top_average = ua_top_average + ua_top_ini;
            va_top_average = va_top_average + va_top_ini;
        end
    end

    ua_skin_all(i) = ua_top_average;
    va_skin_all(i) = va_top_average;

    time_all(i) = 6 * (i - 1);
end

for j = 1:step

% absolute values of relative velocities
    omid = sqrt((uo_cst_all(j) - U) ^ 2 + (vo_cst_all(j) - V) ^ 2);
    amid = sqrt((ua_all(j) - U) ^ 2 + (va_all(j) - V) ^ 2);
    omid_skin = sqrt((uo_skin_all(j) - U) ^ 2 + (vo_skin_all(j) - V) ^ 2);
    amid_skin = sqrt((ua_skin_all(j) - U) ^ 2 + (va_skin_all(j) - V) ^ 2);

% Force due to Air
    Fa_u = rho_air * 0.5 * Ca * dA_a * amid * (ua_all(j) - U) + rho_air * Cda_skin * Ad * amid_skin * (ua_all(j) - U);
    Fa_v = rho_air * 0.5 * Ca * dA_a * amid * (va_all(j) - V) + rho_air * Cda_skin * Ad * amid_skin * (va_all(j) - V);

% Force due to the Ocean
    Fo_u = rho_h2o * 0.5 * Co * dA_o * omid * (uo_cst_all(j) - U) + rho_h2o * Cdo_skin * Ad * omid_skin * (uo_cst_all(j) - U);
    Fo_v = rho_h2o * 0.5 * Co * dA_o * omid * (vo_cst_all(j) - V) + rho_h2o * Cdo_skin * Ad * omid_skin * (vo_cst_all(j) - V);
    
% Coriolis Force & Pressure Gradient Force
    Fcp_u = -(- M * f * (V - vo_cst_all(j)));
    Fcp_v = - M * f * (U - uo_cst_all(j));

% calculate acceleration
    au = (Fcp_u + Fa_u + Fo_u) / M;
    av = (Fcp_v + Fa_v + Fo_v) / M;  
                         
% calculate velocity
  s_u = U * dt + 0.5 * au * dt ^ 2;
  s_v = V * dt + 0.5 * av * dt ^ 2;
  x_all(j) = x;
  y_all(j) = y;
  xx_all(j) = xx;
  yy_all(j) = yy;
  U_all(j) = U;
  V_all(j) = V;
  vel = sqrt(U ^ 2 + V ^ 2);
  vel_all(j) = vel;
  Fa_all(j) = sqrt(Fa_u ^ 2 + Fa_v ^ 2);
  Fo_all(j) = sqrt(Fo_u ^ 2 + Fo_v ^ 2);
  Fcp_all(j) = sqrt(Fcp_u ^ 2 + Fcp_v ^ 2);
    
% new velocity & location
  U = U + au * dt;
  V = V + av * dt;
  xx = xx + s_u;
  yy = yy + s_v;
  if xx > x_rho(x,y)
     x_all(i) = x + 1;
  else
      x_all(i) = x;
  end
  if yy > y_rho(x,y)
     y_all(i) = y + 1;
  else
      y_all(i) = y;
  end

end

for t = 1:step
    if xx_all(t)>=0 && xx_all(t)<=X_lim && yy_all(t)>=0 && yy_all(t)<=Y_lim 
       if mask_zice(x_all(t),y_all(t)) == 0 && mask_land(x_all(t),y_all(t)) == 1
          if depth_icb_under > abs(depths_rho(index_x,index_y,index_z))
             break
          end
       else
           break
       end
    else
        break
    end
end

if depth_icb_under > abs(depths_rho(index_x,index_y,index_z))
    disp(['warning! ground when time = ',num2str(time_all(t-1)),' h']);
    disp(['x = ',num2str(xx_all(t-1)),' y = ',num2str(yy_all(t-1))]);
else
    disp(['warning! stop moving when time = ',num2str(time_all(t-1)),' h']);
    disp(['x = ',num2str(xx_all(t-1)),' y = ',num2str(yy_all(t-1))]);
end

subplot(2,2,1);
plot(time_all(1:t-1),vel_all(1:t-1),'*-');
title('icb velocity - time');
xlabel('time(h)');
ylabel('velocity(m/s)');

subplot(2,2,2);
plot(time_all(1:t-1),U_all(1:t-1),'-r.');
hold on;
plot(time_all(1:t-1),V_all(1:t-1),'-b*');
legend('U','V');
legend('boxoff');
title('icb U&V - time');
xlabel('time(h)');
ylabel('velocity(m/s)');

subplot(2,2,3);
plot(xx_all(1:t-1),yy_all(1:t-1),'*-');
grid;
hold on;
plot(xx_all(1),yy_all(1),'ro');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
% boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,length,length,depth);

subplot(2,2,4);
plot(time_all(1:t-1),Fa_all(1:t-1),'-r.');
hold on;
plot(time_all(1:t-1),Fo_all(1:t-1),'-b*');
hold on;
plot(time_all(1:t-1),Fcp_all(1:t-1),'-y+');
legend('F_ atomas','F_ ocean','F_ pressure&coriolis');
legend('boxoff');
title('F - time');
xlabel('time(h)');
ylabel('force(N)');

% check if within the domain

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
dx_rho = 1./ncread(gname,'pm'); % 96 * 120
dy_rho = 1./ncread(gname,'pn');
h = ncread(gname,'h'); % bathymetry at rho_points (96 * 120)
zice = ncread(gname,'zice'); % ice draft at rho_points (96 * 120)
depths_rho_read = calc_depths_rho(gname); % depths of each level (96 * 120 * 31)
uo_read = ncread(uname,'u'); % meter second-1 (96x119x31x1456)
vo_read = ncread(vname,'v');
sustr_read = ncread(uname,'sustr'); % newton meter-2 (96x119x1456) 
svstr_read = ncread(vname,'svstr');
[X,Y] = meshgrid([1:size(lat_rho,2)],[1:size(lon_rho,1)]);

% indicate the initial location of the particle
x_ini = 5;
y_ini = 10;
x = x_ini;
y = y_ini;

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
f = 2 * w_r * sin(-67*pi/180);
    
% size of the iceberg 
rho_icb = 916.7; % kg m-3
depth = 10; % m
depth_icb_under = depth * (rho_icb / rho_h2o);
depth_icb_above = depth - depth_icb_under;
lenth = 10; % assume width is equal to lenth
dA_o = lenth * depth_icb_under; % m^2 ideally per face 
dA_a = lenth * depth_icb_above; % m^2 ideally per face
Ad = lenth ^ 2;
M = rho_icb * lenth ^ 2 * depth; % kg

% data set
step = 2000;
time_all = zeros(1,step);
U_all = zeros(1,step);
V_all = zeros(1,step);
x_all = zeros(1,step);
y_all = zeros(1,step);
Fa_all = zeros(1,step);
Fo_all = zeros(1,step);
Fc_all = zeros(1,step);
Fcp_all = zeros(1,step);
vel_all = zeros(1,step);
z_all = zeros(1,step);
uo_cst_all = zeros(1,step);
vo_cst_all = zeros(1,step);
ua_all = zeros(1,step);
va_all = zeros(1,step);
uo_cst_surf_all = zeros(1,step);
vo_cst_surf_all = zeros(1,step);
depths_rho = zeros(96,120,31);
s_y = zeros(96,120);
s_x = zeros(96,120);
area_cell_ocean_u = zeros(96,120,31);
area_cell_ocean_v = zeros(96,120,31); 
sustr = zeros(96,120);
sustr = zeros(96,120);
area_cell_air_u = zeros(96,120);
area_cell_air_u = zeros(96,120);

% timestep
dt = 6 * 60 * 60; % s

% create a new matrix of ocean velocity
for XI = 1:95
   for ETA = 1:119
       uo_cst_rho(XI,ETA,:,:) = (uo_read(XI,ETA,:,:) + uo_read(XI+1,ETA,:,:)) / 2;
       vo_cst_rho(XI,ETA,:,:) = (vo_read(XI,ETA,:,:) + vo_read(XI+1,ETA,:,:)) / 2;
   end
end

% do a loop
for i = 1:step

% position of iceberg
    for c = 1:31
        depths_rho(x,y,c) = depths_rho_read(x,y,c) + depths_rho(x,y,c);
        if depth_icb_under < depths_rho(x,y,c)
           break
        end
    end

    for a = x:96
        for b = y:120
            s_y(a,b) = dy_rho(a,b) + s_y(a,b);
            if length < s_y(a,b)
               break
            end
        end
        s_x(a,b) = dx_rho(a,b) + s_x(a,b);
        if length < s_y(a,b)
           break
        end
    end

% area percentage under the ocean at u component
    for b2 = y:(b-1)
        for c2 = 1:(c-1)
            area_cell_ocean_u(x,b2,c2) = dy_rho(x,b2) * depths_rho(x,b2,c2);
        end
    end
    area_total_ocean_u = sum(area_cell_ocean_u(x,:,:));

% area percentage under the ocean at v component
    for a2 = x:(a-1)
        for c2 = 1:(c-1)
            area_cell_ocean_v(a2,y,c2) = dy_rho(a2,y) * depths_rho(a2,y,c2);
        end
    end
    area_total_ocean_v = sum(area_cell_ocean_u(:,y,:));

% average ocean velocity
    uo_average = 0; % meter second-1
    vo_average = 0;
    for a3 = x:(a-1)
        for b3 = y:(b-1)
            for c3 = 1:(c-1)
                uo_ini = uo_cst_rho(a3,b3,c3,step) * (area_cell_ocean_u(a3,b3,c3) / area_total_ocean_u);
                vo_ini = vo_cst_rho(a3,b3,c3,step) * (area_cell_ocean_v(a3,b3,c3) / area_total_ocean_v);
                uo_average = uo_average + uo_ini;
                vo_average = vo_average + vo_ini;
            end
        end
    end

    uo_cst_all(i) = uo_average;
    vo_cst_all(i) = vo_average;

% ua = sqrt(sustr(x,y,step) / (rho_air * Cd)); m s-2

% area percentage above the ocean at u component
    for b4 = y:(b-1)
        area_cell_air_u(x,b4) = dy_rho(x,b4) * depth_icb_above;
    end
    area_total_air_u = sum(area_cell_air_u(x,:));

% area percentage above the ocean at v component
    for a4 = x:(a-1)
        area_cell_air_v(a4,y) = dy_rho(a4,y) * depth_icb_above;
    end
    area_total_air_v = sum(area_cell_air_u(:,y));

% average atmospheric velocity
    ua_average = 0; % meter second-1
    va_average = 0;
    for a5 = x:(a-1)
        for b5 = y:(b-1)
            ua_ini = sqrt(sustr_read(x,y,step) / (rho_air * Cd)) * (area_cell_air_u(a5,b5) / area_total_air_u);
            va_ini = sqrt(sustr_read(x,y,step) / (rho_air * Cd)) * (area_cell_air_v(a5,b5) / area_total_air_v);
            ua_average = ua_average + ua_ini;
            va_average = va_average + va_ini;
        end
    end
    ua_all(i) = ua_average;
    va_all(i) = va_average;

    time_all(i) = i / (60 * 60);
end

for j = 1:step

% absolute values of relative velocities
    omid = sqrt((uo_cst_all(j) - U) ^ 2 + (vo_cst_all(j) - V) ^ 2);
    amid = sqrt((ua_all(j) - U) ^ 2 + (va_all(j) - V) ^ 2);
    omid_skin = sqrt((uo_cst_all(j) - U) ^ 2 + (vo_cst_all(j) - V) ^ 2);
    amid_skin = sqrt((ua_all(j) - U) ^ 2 + (va_all(j) - V) ^ 2);

% Force due to Air
    Fa_u = rho_air * 0.5 * Ca * dA_a * amid * (ua_all(j) - U) + rho_air * Cda_skin * Ad * amid_skin * (ua_all(j) - U);
    Fa_v = rho_air * 0.5 * Ca * dA_a * amid * (va_all(j) - V) +  rho_air * Cda_skin * Ad * amid_skin * (va_all(j) - V);

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
  U_all(j) = U;
  V_all(j) = V;
  vel = sqrt(U ^ 2 + V ^ 2);
  vel_all(j) = vel;
  Fa_all(j) = sqrt(Fa_u ^ 2 + Fa_v ^ 2);
  Fo_all(j) = sqrt(Fo_u ^ 2 + Fo_v ^ 2);
  Fc_all(j) = sqrt(Fc_u ^ 2 + Fc_v ^ 2);
  Fcp_all(j) = sqrt(Fcp_u ^ 2 + Fcp_v ^ 2);
    
% new velocity & location
  U = U + au * dt;
  V = V + av * dt;
  x = x + s_u;
  y = y + s_v;

end

for t = 1:step
    if x_all(t)>=0 && x_all(t)<=X_lim && y_all(t)>=0 && y_all(t)<=Y_lim 
       if depth_icb_under > abs(depths_rho(a,b,c))
          break
       end
    else
        break
    end
end 

if depth_icb_under > abs(depths_rho(a,b,c))
    disp(['warning! ground when time = ',num2str(time_all(t-1)),' h']);
    disp(['x = ',num2str(x_all(t-1)),' y = ',num2str(y_all(t-1))]);
else
    disp(['warning! out of domain when time = ',num2str(time_all(t-1)),' h']);
    disp(['x = ',num2str(x_all(t-1)),' y = ',num2str(y_all(t-1))]);
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
mesh(X,Y,Z);
hold on;
zlim([-50,0]);
plot(x_all(1:t-1),y_all(1:t-1),'*-');
grid;
hold on;
plot(x_all(1),y_all(1),'ro');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
zlabel('depth(m)');
boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,lenth,width,depth);

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

%{
figure, mesh(X,Y,Z);
hold on;
zlim([-50,5]);
plot(x_all(1:t-1),y_all(1:t-1),'*-');
grid;
hold on;
plot(x_all(1),y_all(1),'ro');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
zlabel('depth(m)');
boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,lenth,width,depth);
%}

% check if within the domain

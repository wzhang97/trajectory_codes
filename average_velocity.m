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
strname = '/ds/projects/mertz/mdl/misom020_2/mer_his_0001.nc';
lat_rho = ncread(gname,'lat_rho'); % 96 * 120
lon_rho = ncread(gname,'lon_rho');
x_rho = ncread(gname,'x_rho',[1 1],[94 118]);
y_rho = ncread(gname,'y_rho',[1 1],[94 118]);
[dx_rho,dy_rho,dz_rho] = calc_rhogrid_dxdydz(gname); % calcualte the length of each edges of a cell
depths_rho = calc_depths_rho(gname); % depths of each level (96 * 120 * 31)
bathy = ncread(gname,'h');
mask_zice = ncread(gname,'mask_zice'); % mask on ice-shelf rho-points (option(1) = 'land', option(0) = 'water')
mask_land = ncread(gname,'mask_rho'); % mask on RHO-points (option(0) = 'land', option(1) = 'water')
uo_read = ncread(uname,'u'); % meter second-1 (95x119x31x1456)
vo_read = ncread(vname,'v');
sustr_read = ncread(uname,'sustr'); % newton meter-2 (95x120x1456) 
svstr_read = ncread(vname,'svstr'); % (96x119x1456)
Cs_r = ncread(strname,'Cs_r');
%[u_rho,v_rho] = velocity2rho(gname,depths_rho,uo_read,vo_read,0); % calculate velocities on rho-grid
load /v_mertz/iceberg_project/wzhang/trajectory_codes/vel_rho.mat;
[X,Y] = meshgrid([1:size(x_rho,2)],[1:size(y_rho,1)]);
uo_cst_rho = u_rho;
vo_cst_rho = v_rho;
end_layer = length(Cs_r);

% indicate the initial location of the particle
x_ini = 40;
y_ini = 60;
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
rho_h2o = 1027.5; % seawater density kg m-3

% atmospheric coeff
rho_air = 1.225; % air density kg m-3
Ca = 0.4; % dimensionless coefficient of resistance
Cda_skin = 0.0022; % air-iceberg friction coeff
Cd = 1.25e-3;  % dimensionless air-sea friction coeff

% coriolis coeff
w_r = 7.2921e-5;
f = 2 * w_r * sin(-67 * pi / 180);
    
% size of the iceberg 
rho_icb = 850; % kg m-3
depth = 100; % m
depth_icb_under = depth * (rho_icb / rho_h2o);
depth_icb_above = depth - depth_icb_under;
len = 4000; % assume width is equal to length
dA_o = len * depth_icb_under; % m^2 ideally per face 
dA_a = len * depth_icb_above; % m^2 ideally per face
Ad = len ^ 2;
M = rho_icb * len ^ 2 * depth; % kg

% data set
step = 5;
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
uo_cst_all = zeros(1,step);
vo_cst_all = zeros(1,step);
ua_all = zeros(1,step);
va_all = zeros(1,step);
au_all = zeros(1,step);
av_all = zeros(1,step);
uo_skin_all = zeros(1,step);
vo_skin_all = zeros(1,step);
ua_skin_all = zeros(1,step);
va_skin_all = zeros(1,step);
Fa_drag_uall = zeros(1,step);
Fa_skin_uall = zeros(1,step);
Fa_drag_vall = zeros(1,step);
Fa_skin_vall = zeros(1,step);
Fo_drag_uall = zeros(1,step);
Fo_skin_uall = zeros(1,step);
Fo_drag_vall = zeros(1,step);
Fo_skin_vall = zeros(1,step);
Fa_u_all = zeros(1,step);
Fa_v_all = zeros(1,step);
Fo_u_all = zeros(1,step);
Fo_v_all = zeros(1,step);
Fcp_u_all = zeros(1,step);
Fcp_v_all = zeros(1,step);



% timestep
dt = 6 * 60 * 60; % s

% create a new matrix of wind velocity 
% ua = sqrt(sustr(x,y,step) / (rho_air * Cd)); m s-2
%{
for XI = 1:94 
    for ETA = 1:118
        sustr_rho(XI,ETA,:) = 0.5 .* (sustr_read(XI,ETA,:) + sustr_read(XI+1,ETA,:));
        u_wind = sign(sustr_rho) .* sqrt(abs(sustr_rho) / (rho_air * Cd));
    end
end

for XI = 1:94
    for ETA = 1:118
        svstr_rho(XI,ETA,:) = 0.5 .* (svstr_read(XI,ETA,:) + svstr_read(XI,ETA+1,:));
        v_wind = sign(svstr_rho) .* sqrt(abs(svstr_rho) / (rho_air * Cd));
    end
end
%}
% do a loop
for i = 1:step

% position of iceberg
    tmp_z = depths_rho(x,y,:);
    [closest_diff_z, index_z] = min(abs((-tmp_z) - depth_icb_under)); % closest_diff is the difference between the bottom and the closest layer, index is the number of the closest layer
    len_rho_x = 0;
    len_rho_y = 0;

    for a = x:94        
        len_rho_x = dx_rho(a,y) + len_rho_x;
        if len < len_rho_x 
           if  len > (len_rho_x - 0.5 * dx_rho(a,y))
               index_x = a + 1;
               break
           else 
               index_x = a;
               break
           end
        end
    end
           
    for b = y:118
        len_rho_y = dy_rho(x,b) + len_rho_y;
        if len < len_rho_y
           if  len > (len_rho_y - 0.5 * dy_rho(x,b))
               index_y = b + 1;
               break
           else
               index_y = b;
               break
           end
        end
    end
% area percentage under the ocean at u component(side)
    Num_cell_hori_u = index_y - y; 
    Num_cell_vert = 31 - index_z + 1;
    dy_rho_repM = repmat(dy_rho,[1,1,31]);
    area_cell_ocean_u = dy_rho_repM(x,y:index_y,index_z:end_layer) .* dz_rho(x,y:index_y,index_z:end_layer);
    area_cell_ocean_u = squeeze(area_cell_ocean_u); % to remove the ﬁrst dimension, check the evolution of the size of the matrix when you apply that 
    area_total_ocean_u = sum(sum(area_cell_ocean_u));

% area percentage under the ocean at v component(side)
    Num_cell_hori_v = index_x - x;
    dx_rho_repM = repmat(dx_rho,[1,1,31]);
    area_cell_ocean_v = dx_rho_repM(x:index_x,y,index_z:end_layer) .* dz_rho(x:index_x,y,index_z:end_layer);
    area_cell_ocean_v = squeeze(area_cell_ocean_v);
    area_total_ocean_v = sum(sum(area_cell_ocean_v));

% average ocean velocity on side
    uo_weighted = squeeze(uo_cst_rho(x,y:index_y,index_z:end_layer,i)) .* (area_cell_ocean_u(:,:) ./ area_total_ocean_u);
    uo_average = sum(uo_weighted(:));

    vo_weighted = squeeze(vo_cst_rho(x:index_x,y,index_z:end_layer,i)) .* (area_cell_ocean_v(:,:) ./ area_total_ocean_v);
    vo_average = sum(vo_weighted(:));

    uo_cst_all(i) = uo_average;
    vo_cst_all(i) = vo_average;

% area percentage at the bottom
    area_cell_bottom = dx_rho(x:index_x,y:index_y) .* dy_rho(x:index_x,y:index_y);
    area_total_bottom = sum(sum(area_cell_bottom));

% average ocean velocity at bottom
    uo_bottom_weighted = squeeze(uo_cst_rho(x:index_x,y:index_y,index_z,i)) .* (area_cell_bottom(:,:) ./ area_total_bottom);
    uo_bottom_average = sum(uo_bottom_weighted(:));

    vo_bottom_weighted = squeeze(vo_cst_rho(x:index_x,y:index_y,index_z,i)) .* (area_cell_bottom(:,:) ./ area_total_bottom);
    vo_bottom_average = sum(vo_bottom_weighted(:));

    uo_skin_all(i) = uo_bottom_average;
    vo_skin_all(i) = vo_bottom_average;

% area percentage above the ocean at u component
    area_cell_air_u = dy_rho(x,y:index_y) .* depth_icb_above;

    area_total_air_u = sum(sum(area_cell_air_u));

% area percentage above the ocean at v component
    area_cell_air_v = dx_rho(x:index_x,y) .* depth_icb_above;

    area_total_air_v = sum(sum(area_cell_air_v));

% average atmospheric velocity on side
    ua_weighted = squeeze(u_wind(x,y:index_y,i)) .* (area_cell_air_u(:,:) ./ area_total_air_u);
    ua_average = sum(ua_weighted(:));

    va_weighted = squeeze(v_wind(x:index_x,y,i)) .* (area_cell_air_v(:,:) ./ area_total_air_v);
    va_average = sum(va_weighted(:));

    ua_all(i) = ua_average;
    va_all(i) = va_average;

% average ocean velocity at top
    ua_top_weighted = squeeze(u_wind(x:index_x,y:index_y,i)) .* (area_cell_bottom(:,:) ./ area_total_bottom);
    ua_top_average = sum(ua_top_weighted(:));

    va_top_weighted = squeeze(v_wind(x:index_x,y:index_y,i)) .* (area_cell_bottom(:,:) ./ area_total_bottom);
    va_top_average = sum(va_top_weighted(:));

    ua_skin_all(i) = ua_top_average;
    va_skin_all(i) = va_top_average;

    time_all(i) = 6 * (i - 1);

% absolute values of relative velocities
    omib = sqrt((uo_cst_all(i) - U) ^ 2 + (vo_cst_all(i) - V) ^ 2);
    amib = sqrt((ua_all(i) - U) ^ 2 + (va_all(i) - V) ^ 2);
    omib_skin = sqrt((uo_skin_all(i) - U) ^ 2 + (vo_skin_all(i) - V) ^ 2);
    amib_skin = sqrt((ua_skin_all(i) - U) ^ 2 + (va_skin_all(i) - V) ^ 2);

% Force due to Air
    Fa_drag_u = rho_air * 0.5 * Ca * dA_a * amib * (ua_all(i) - U);
    Fa_skin_u = rho_air * Cda_skin * Ad * amib * (ua_all(i) - U);
    Fa_drag_v = rho_air * 0.5 * Ca * dA_a * amib * (va_all(i) - V);
    Fa_skin_v = rho_air * Cda_skin * Ad * amib * (va_all(i) - V);

    Fa_u = rho_air * 0.5 * Ca * dA_a * amib * (ua_all(i) - U) + rho_air * Cda_skin * Ad * amib * (ua_all(i) - U);
    Fa_v = rho_air * 0.5 * Ca * dA_a * amib * (va_all(i) - V) + rho_air * Cda_skin * Ad * amib * (va_all(i) - V);

% Force due to the Ocean
%    Fo_u = rho_h2o * 0.5 * Co * dA_o * omib * (uo_cst_all(i) - U) + rho_h2o * Cdo_skin * Ad * omib_skin * (uo_cst_all(i) - U);
%    Fo_v = rho_h2o * 0.5 * Co * dA_o * omib * (vo_cst_all(i) - V) + rho_h2o * Cdo_skin * Ad * omib_skin * (vo_cst_all(i) - V);
    Fo_drag_u = rho_h2o * 0.5 * Co * dA_o * omib * uo_cst_all(i);
    Fo_skin_u = rho_h2o * Cdo_skin * Ad * omib_skin * uo_skin_all(i);
    Fo_drag_v = rho_h2o * 0.5 * Co * dA_o * omib * vo_cst_all(i);
    Fo_skin_v = rho_h2o * Cdo_skin * Ad * omib_skin * vo_skin_all(i);

    Fo_u = rho_h2o * 0.5 * Co * dA_o * omib * uo_cst_all(i) + rho_h2o * Cdo_skin * Ad * omib_skin * uo_skin_all(i);
    Fo_v = rho_h2o * 0.5 * Co * dA_o * omib * vo_cst_all(i) + rho_h2o * Cdo_skin * Ad * omib_skin * vo_skin_all(i);

% Coriolis Force & Pressure Gradient Force
    Fcp_u = -(- M * f * (V - vo_cst_all(i)));
    Fcp_v = - M * f * (U - uo_cst_all(i));

% calculate acceleration
    au = (Fcp_u + Fa_u + Fo_u) / M;
    av = (Fcp_v + Fa_v + Fo_v) / M; 
    au_all(i) = au;
    av_all(i) = av;
% calculate velocity
    s_u = U * dt + 0.5 * au * dt ^ 2;
    s_v = V * dt + 0.5 * av * dt ^ 2;
    xx_all(i) = xx;
    yy_all(i) = yy;
    x_all(i) = x;
    y_all(i) = y;
    U_all(i) = U;
    V_all(i) = V;
    vel = sqrt(U ^ 2 + V ^ 2);
    vel_all(i) = vel;
    Fa_drag_uall(i) = Fa_drag_u;
    Fa_skin_uall(i) = Fa_skin_u;
    Fa_drag_vall(i) = Fa_drag_v;
    Fa_skin_vall(i) = Fa_skin_v;
    Fo_drag_uall(i) = Fo_drag_u;
    Fo_skin_uall(i) = Fo_skin_u;
    Fo_drag_vall(i) = Fo_drag_v;
    Fo_skin_vall(i) = Fo_skin_v;
    Fa_u_all(i) = Fa_u;
    Fa_v_all(i) = Fa_v;
    Fo_u_all(i) = Fo_u;
    Fo_v_all(i) = Fo_v;
    Fcp_u_all(i) = Fcp_u;
    Fcp_v_all(i) = Fcp_v;
    Fa_all(i) = sqrt(Fa_u ^ 2 + Fa_v ^ 2);
    Fo_all(i) = sqrt(Fo_u ^ 2 + Fo_v ^ 2);
    Fcp_all(i) = sqrt(Fcp_u ^ 2 + Fcp_v ^ 2);
    
% new velocity & location
    U = U + au * dt;
    V = V + av * dt;
    xx = xx + s_u;
    yy = yy + s_v;
    y = griddata(x_rho, y_rho, X, double(xx), double(yy), 'nearest');
    x = griddata(x_rho, y_rho, Y, double(xx), double(yy), 'nearest');

    if (x > size(x_rho,1)) || (y > size(y_rho,2))    
       disp(['outside of the domain! when time = ',num2str(time_all(i-1)),' h']);
       disp(['x = ',num2str(xx_all(i-1)),' y = ',num2str(yy_all(i-1))]);
       break
    elseif mask_zice(x,y) == 1
       disp(['iceberg on the ice! when time = ',num2str(time_all(i-1)),' h']);
       disp(['x = ',num2str(xx_all(i-1)),' y = ',num2str(yy_all(i-1))]);
       break 
    elseif mask_land(x,y) == 0    
       disp(['iceberg on land! when time = ',num2str(time_all(i-1)),' h']);
       disp(['x = ',num2str(xx_all(i-1)),' y = ',num2str(yy_all(i-1))]);
       break 
    elseif abs(bathy(x,y)) < abs(depth_icb_under)
       disp(['grounded! when time = ',num2str(time_all(i-1)),' h']);
       disp(['x = ',num2str(xx_all(i-1)),' y = ',num2str(yy_all(i-1))]);
       break 
    end

end

%{
subplot(2,2,1);
plot(time_all(1:i-1),vel_all(1:i-1),'*-');
title('icb velocity - time');
xlabel('time(h)');
ylabel('velocity(m/s)');

subplot(2,2,2);
plot(time_all(1:i-1),U_all(1:i-1),'-r.');
hold on;
plot(time_all(1:i-1),V_all(1:i-1),'-b*');
legend('U','V');
legend('boxoff');
title('icb U&V - time');
xlabel('time(h)');
ylabel('velocity(m/s)');

subplot(2,2,3);
plot(xx_all(1:t-1),yy_all(1:i-1),'*-');
grid;
hold on;
plot(xx_all(1),yy_all(1),'ro');
plot(xx_all(i-1),yy_all(i-1),'bo');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
% boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,len,len,depth);

subplot(2,2,4);
plot(time_all(1:i-1),Fa_all(1:i-1),'-r.');
hold on;
plot(time_all(1:i-1),Fo_all(1:i-1),'-b*');
hold on;
plot(time_all(1:i-1),Fcp_all(1:i-1),'-y+');
legend('F_ atomas','F_ ocean','F_ pressure&coriolis');
legend('boxoff');
title('F - time');
xlabel('time(h)');
ylabel('force(N)');
%}
% check if within the domain

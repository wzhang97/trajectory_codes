%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Grounded model
%
%    Author: weiyan ZHANG
%    Created: Oct 2018
%    Updated:
%    Last Modif by:
%    Last Modif on: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% load grid information from the model
%{
gname = 'C:\Users\61414\Downloads\misom020_grd.nc';
lat_rho = ncread(gname,'lat_rho');
lon_rho = ncread(gname,'lon_rho');
dx = ones(size(lat_rho,1),size(lat_rho,2));
dy = ones(size(lat_rho,1),size(lat_rho,2));
[xx yy] = meshgrid([1:size(lat_rho,2)],[1:size(lon_rho,1)]); 
%}

%{
% calculate dx, dy from model in m
pm = ncread(gname,'pm');
pn = ncread(gname,'pn');
dx = 1./pm;
dy = 1./pn;
%}

% Creates 2d grid of same size as model
X_lim = 10000;
Y_lim = 10000;
[X,Y] = meshgrid([0:X_lim],[0:Y_lim]);
Z = 10 - X ;

% indicate the initial location of the particle
xx_ini = 350;
yy_ini = 500;

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

% timestep
dt = 60; % s

% do a loop    
xx_ = xx_ini;
yy_ = yy_ini;

for i=1:step
    
    theta = (i/step) * 8 * pi;
% ocean velocity
    uo_cst = 0.1 * sin(theta); % m/s
    vo_cst = 0.15 * sin(theta);

% atmospheric velocity (at 10m)
% the wind in the model will be given as wind stress -- do the conversion
% tau_wind = rho_air*Cd*U^2 (tau_wind is sustr or svstr in model)

    sustr = abs(0.05 * sin(theta)); %newton meter-2
    svstr = abs(0.035 * sin(theta));
    ua = sqrt(sustr / (rho_air * Cd)); % m s-2
    va = sqrt(svstr / (rho_air * Cd)); % m s-2
    
    time_all(i) = i/60;
    uo_cst_all(i) = uo_cst;
    vo_cst_all(i) = vo_cst;
    ua_all(i) = ua;
    va_all(i) = va;
end

for j = 1:step

% absolute values of relative velocities
    omid = sqrt((uo_cst_all(j) - U) ^ 2 + (vo_cst_all(j) - V) ^ 2);
    amid = sqrt((ua_all(j) - U) ^ 2 + (va_all(j) - V) ^ 2);
    omid_skin = sqrt((uo_cst_all(j) - U) ^ 2 + (vo_cst_all(j) - V) ^ 2);
    amid_skin = sqrt((ua_all(j) - U) ^ 2 + (va_all(j) - V) ^ 2);

% Force due to the Atomasphere
    Fa_u = rho_air * 0.5 * Ca * dA_a * amid * (ua_all(j) - U) + rho_air * Cda_skin * Ad * amid_skin * (ua_all(j) - U);
    Fa_v = rho_air * 0.5 * Ca * dA_a * amid * (va_all(j) - V) +  rho_air * Cda_skin * Ad * amid_skin * (va_all(j) - V);

% Force due to the Ocean
    Fo_u = rho_h2o * 0.5 * Co * dA_o * omid * (uo_cst_all(j) - U) + rho_h2o * Cdo_skin * Ad * omid_skin * (uo_cst_all(j) - U);
    Fo_v = rho_h2o * 0.5 * Co * dA_o * omid * (vo_cst_all(j) - V) + rho_h2o * Cdo_skin * Ad * omid_skin * (vo_cst_all(j) - V);

% Coriolis Force
    Fc_u = - (M * f * V);
    Fc_v = - M * f * U;
    
% Coriolis Force & Pressure Gradient Force
    Fcp_u = -(- M * f * (V - vo_cst_all(j)));
    Fcp_v = - M * f * (U - uo_cst_all(j));
%{    
% calculate acceleration
    au = (Fc_u + Fa_u + Fo_u) / M;
    av = (Fc_v + Fa_v + Fo_v) / M;  
%}
% calculate acceleration
    au = (Fcp_u + Fa_u + Fo_u) / M;
    av = (Fcp_v + Fa_v + Fo_v) / M;  
    
% calculate velocity
  s_u = U * dt + 0.5 * au * dt ^ 2;
  s_v = V * dt + 0.5 * av * dt ^ 2;
  x_all(j) = xx_;
  y_all(j) = yy_;
  U_all(j) = U;
  V_all(j) = V;
  vel = sqrt(U ^ 2 + V ^ 2);
  vel_all(j) = vel;
  Fa_all(j) = sqrt(Fa_u ^ 2 + Fa_v ^ 2);
  Fo_all(j) = sqrt(Fo_u ^ 2 + Fo_v ^ 2);
  Fc_all(j) = sqrt(Fc_u ^ 2 + Fc_v ^ 2);
  Fcp_all(j) = sqrt(Fcp_u ^ 2 + Fcp_v ^ 2);
  
  z = 10 + xx_;
  z_all(j) = z;
  
% new velocity & location
  U = U + au * dt;
  V = V + av * dt;
  xx_ = xx_ + s_u;
  yy_ = yy_ + s_v;
end

for t = 1:step
    if x_all(t)>=0 && x_all(t)<=X_lim && y_all(t)>=0 && y_all(t)<=Y_lim 
        if depth_icb_under > abs(z_all(t))
            break
        end
    else
        break
    end
end 

if abs(depth_icb_under) > abs(z)
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
title('icb U & V - time');
xlabel('time(h)');
ylabel('velocity(m/s)');

subplot(2,2,3);
%mesh(X,Y,Z);
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
%boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,lenth,width,depth);

subplot(2,2,4);
plot(time_all(1:t-1),Fa_all(1:t-1),'-r.');
hold on;
plot(time_all(1:t-1),Fo_all(1:t-1),'-b*');
hold on;
plot(time_all(1:t-1),Fc_all(1:t-1),'-y+');
hold on;
plot(time_all(1:t-1),Fcp_all(1:t-1),'-g+');
legend('F_ atomas','F_ ocean','F_ coriolis','F_ pressure&coriolis');
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
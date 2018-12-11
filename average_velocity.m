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

% load grid information from the model
gname = 'C:\Users\61414\Downloads\misom020_grd.nc';
lat_rho = ncread(gname,'lat_rho');
lon_rho = ncread(gname,'lon_rho');
%{
% calculate dx, dy from model in m
pm = ncread(gname,'pm');
pn = ncread(gname,'pn');
dx = 1./pm;
dy = 1./pn;
%}
dx = ones(size(lat_rho,1),size(lat_rho,2));
dy = ones(size(lat_rho,1),size(lat_rho,2));

% Creates 2d grid of same size as model
X_lim = 1000;
Y_lim = 1000;
[xx yy] = meshgrid([1:size(lat_rho,2)],[1:size(lon_rho,1)]); 
[X,Y] = meshgrid([0:X_lim],[0:Y_lim]);
Z = -300 + X ;
%Z = - (xx_,yy_);

% indicate the initial location of the particle
xx_ini = 50;
yy_ini = 150;
xx_ = xx_ini;
yy_ = yy_ini;

% set variable useful for the trajectory equation
U = 0;
V = 0;

% size of the iceberg
depth = 10; % m
depth_icb_under = 9;
depth_icb_above = 1;
lenth = 10;
width = 10;
dA_o = width * depth_icb_under; % m^2 ideally per face 
dA_a = width * depth_icb_above; % m^2 ideally per face
Ad = lenth * width;
M = 916.7 * lenth * width * depth; % kg
                                   % density is 916.7 kg m-3
             
step = 500;
time_all = zeros(1,step);
U_all = zeros(1,step);
V_all = zeros(1,step);
x_all = zeros(1,step);
y_all = zeros(1,step);
Fa_all = zeros(1,step);
Fo_all = zeros(1,step);
Fc_all = zeros(1,step);
vel_all = zeros(1,step);
z_all = zeros(1,step);
uo_cst_all = zeros(1,step);
vo_cst_all = zeros(1,step);
ua_all = zeros(1,step);
va_all = zeros(1,step);

% oceanic coeff
Co = 0.85; % dimensionless coefficient of resistance
Cdo_skin = 0.0055; % ocean-iceberg friction coeff
rho_h2o = 1028; % seawater density kg m-3

% atmospheric coeff
Ca = 0.4; % dimensionless coefficient of resistance
Cda_skin = 0.0022; % air-iceberg friction coeff
rho_air = 1.225; % air density kg m-3

% atmospheric velocity (at 10m)
% the wind in the model will be given as wind stress -- do the conversion
% tau_wind = rho_air*Cd*U^2 (tau_wind is sustr or svstr in model)
sustr = 0.05; %newton meter-2
svstr = 0.05;
Cd = 1.25e-3;  % dimensionless air-sea friction coeff
ua= sqrt(sustr / (rho_air * Cd)); % m s-2
va= sqrt(svstr / (rho_air * Cd)); % m s-2

% timestep
dt = 60; % s

% coriolis coeff
w_r = 7.2921e-5;
f = - 2 * w_r * sin(67*pi/180);

% do a loop    
for i=1:step
    uv=sqrt(U^2+V^2);
    
% coastal current
    uo_cst_surf = 0.0015 * 0.1 * yy_;
    vo_cst_surf = 0.0015 * 0.1 * yy_;

% ocean velocity
    uo_cst_ans = 0;
    vo_cst_ans = 0;
    level = depth_icb_under / 1;
for depth = 0:level
    uo_cst_ini = uo_cst_surf - 0.005 * depth; % m/s
    vo_cst_ini = vo_cst_surf - 0.005 * depth;

    uo_cst_ans = uo_cst_ans + uo_cst_ini;
    vo_cst_ans = vo_cst_ans + vo_cst_ini;
end
% average ocean velocity
    uo_cst = uo_cst_ans / level;
    vo_cst = vo_cst_ans / level;

% absolute values of relative velocities
        omid = sqrt((uo_cst - U) ^ 2 + (vo_cst - V) ^ 2);
        amid = sqrt((ua - U) ^ 2 + (va - V) ^ 2);
        omid_skin = sqrt((uo_cst - U) ^ 2 + (vo_cst - V) ^ 2);
        amid_skin = sqrt((ua - U) ^ 2 + (va - V) ^ 2);

% Force due to Air
        Fa_u = rho_air * 0.5 * Ca * dA_a * amid * (ua - U) + rho_air * Cda_skin * Ad * amid_skin * (ua - U);
        Fa_v = rho_air * 0.5 * Ca * dA_a * amid * (va - U) +  rho_air * Cda_skin * Ad * amid_skin * (va - V);

% Force due to the Ocean
        Fo_u = rho_h2o * 0.5 * Co * dA_o * omid * (uo_cst - U) + rho_h2o * Cdo_skin * Ad * omid_skin * (uo_cst - U);
        Fo_v = rho_h2o * 0.5 * Co * dA_o * omid * (vo_cst - U) + rho_h2o * Cdo_skin * Ad * omid_skin * (vo_cst - V);
    
        Fcoriolis = sqrt((-M * f * (V - vo_cst))^2+(M * f * (U - uo_cst))^2);   
        au = ((-M * f * (V - vo_cst)) + Fa_u + Fo_u) / M;
        av = (-(-M * f * (U - uo_cst)) + Fa_v + Fo_v) / M;
                     
% calculate velocity
  time_all(i) = i/60;
  s_u = U*dt+0.5*au*dt^2;
  s_v = V*dt+0.5*av*dt^2;
  x_all(i) = xx_;
  y_all(i) = yy_;
  U_all(i) = U;
  V_all(i) = V;
  vel = sqrt(U^2+V^2);
  vel_all(i) = vel;
  Fa_all(i) = sqrt(Fa_u ^ 2 + Fa_v ^ 2);
  Fo_all(i) = sqrt(Fo_u ^ 2 + Fo_v ^ 2);
  Fc_all(i) = Fcoriolis;
  uo_cst_all(i) = uo_cst;
  vo_cst_all(i) = vo_cst;
  ua_all(i) = ua;
  va_all(i) = va;
  
  z = -300 + xx_;
  z_all(i) = z;

  if xx_>=0 && xx_<=X_lim && yy_>=0 && yy_<=Y_lim 
      if abs(depth_icb_under) > abs(z)
          break
      end
  else
      break
  end
  
% new velocity
  U = U + au * dt;
  V = V + av * dt; 
  xx_ = xx_+s_u;
  yy_ = yy_+s_v;
end

if abs(depth_icb_under) > abs(z)
    disp(['warning! ground when time = ',num2str(time_all(i-1)),' h']);
    disp(['x = ',num2str(x_all(i-1)),' y = ',num2str(y_all(i-1))]);
else
    disp(['warning! out of domain when time = ',num2str(time_all(i-1)),' h']);
    disp(['x = ',num2str(x_all(i-1)),' y = ',num2str(y_all(i-1))]);
end

subplot(2,2,1);
plot(vel_all(1:i-1),'*-');
title('velocity - time');
xlabel('time(min)');
ylabel('velocity(m/s)');

subplot(2,2,2);
plot(U_all(1:i-1),'-r.');
hold on;
plot(V_all(1:i-1),'-b*');
legend('U','V');
legend('boxoff');
title('U&V - time');
xlabel('time(min)');
ylabel('velocity(m/s)');

subplot(2,2,3);
mesh(X,Y,Z);
hold on;
zlim([-50,0]);
plot(x_all(1:i),y_all(1:i),'*-');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
zlabel('depth(m)');

subplot(2,2,4);
plot(Fa_all(1:i-1),'-r.');
hold on;
plot(Fo_all(1:i-1),'-b*');
hold on;
plot(Fc_all(1:i-1),'-y+');
legend('F_ atomas','F_ ocean','F_ coriolis');
legend('boxoff');
title('F - time');
xlabel('time(min)');
ylabel('force(N)');

% check if within the domain

% keep moving in time

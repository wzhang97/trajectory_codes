load C:\Users\61414\Desktop\average_velocity.mat;

x1 = [1:size(bathy,1)];
x2 = repmat(x1,[1,size(bathy,2)]);
y1 = [1:size(bathy,2)];
y2 = kron(y1,ones(size(bathy,1),1));
y3 = y2(:)';
z1 = - bathy(:);
x_rho1 = x_rho(1,:);
y_rho1 = y_rho(:,1);
[X,Y] = meshgrid([1:size(bathy,1)],[1:size(bathy,2)]);
Z=griddata(x2,y3,z1,X,Y,'v4');
surf(X,Y,Z);
shading interp;
colormap(jet);
colorbar;
if mask_zice == 1
    pathc(mask_land, 'g');
end
xlabel('xx(m)');
ylabel('yy(m)');
set(gca,'xtick', x_rho1);
set(gca,'ytick', y_rho1);

plot(xx_all(1:i-1),yy_all(1:i-1),'*-');
hold on;
plot(xx_all(1),yy_all(1),'ro');
hold on;
plot(xx_all(i-1),yy_all(i-1),'bo');
title('trajectory of iceberg');
xlabel('xx(m)');
ylabel('yy(m)');
% boxplot3(x_all(t)-5,y_all(t)-5,-depth_icb_under,len,len,depth);

clear; clc; close all;

load("Lidar_cleaned_1M_pnts.mat");
lidar_pt = pointCloud(xyz_P'); %convert the matrix data into a point cloud, data needs to be Mx3


pcwrite(lidar_pt, 'Lidar.pcd') %write the point cloud in a usable format for the pcread


pt_lidar = pcread('Lidar.pcd'); %display the point cloud data

figure;
pcshow(pt_lidar)
title("PC lidar data");
xlabel("X");
ylabel("Y");
zlabel("Z");
xlim([-2 2]);
ylim([-2 2]);
zlim([0 4]);

figure;
plot3(xyz_P(1,:), xyz_P(2,:), xyz_P(3,:), '.'); %display data as a matrix;
title("Plot3 lidar data");
xlabel("X");
ylabel("Y");
zlabel("Z");
xlim([-2 2]);
ylim([-2 2]);
zlim([0 4]);


%setup for pcfitplane
max_distance = 0.02; %maximum point to plane distance (2 cm) for plane fitting
ref_vector = [0,0,1]; %normal vector to the plane
max_angular_dist = 5; %maximum angular distance

%Detect the first plane 
[model, inlierIndicies, outlierIndicies] = pcfitplane(pt_lidar, max_distance, ref_vector, max_angular_dist);
plane = select(pt_lidar, inlierIndicies);
rest_pt_lidar = select(pt_lidar, outlierIndicies);

figure;
subplot(2, 1, 1);
pcshow(plane);
title("Cieling");
xlabel("X");
ylabel("Y");
zlabel("Z");

%compute the residual of the inliers to the model
d = abs(xyz_P(:, inlierIndicies)' - plane.Location);
subplot(2, 1, 2);
histogram(d);
title("Residual Histogram", 'Color', 'white')
xlabel("bin", 'Color', 'white');
ylabel("frequency", 'Color', 'white');
ax = gca;
ax.YColor = 'w';
ax.XColor = 'w';


%compute and display the outliers

ref_vector = [0,1,0]; %normal vector to the plane
[model, inlierIndicies, outlierIndicies] = pcfitplane(pt_lidar, max_distance, ref_vector, max_angular_dist);
plane = select(pt_lidar, inlierIndicies);
rest_pt_lidar = select(pt_lidar, outlierIndicies);
figure;
pcshow(plane);
title("outliers");



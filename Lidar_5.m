clear; clc; close all;

load("Lidar_cleaned_1M_pnts.mat");
lidar_pt = pointCloud(xyz_P'); %convert the matrix data into a point cloud, data needs to be Mx3


pcwrite(lidar_pt, 'Lidar.pcd') %write the point cloud in a usable format for the pcread


pt_lidar = pcread('Lidar.pcd'); %display the point cloud data

%setup for pcfitplane
max_distance = 0.02; %maximum point to plane distance (2 cm) for plane fitting
ref_vector = [0,0,1]; %normal vector to the plane
max_angular_dist = 5; %maximum angular distance


prev_plane_norm = zeros(3, 1); %pre allocate space to collect plane models
i = 1;
while true 
    [model, inlierIndicies, outlierIndicies] = pcfitplane(pt_lidar, max_distance);
    plane = select(pt_lidar, inlierIndicies);

    figure;
    pcshow(plane);
    title("plane: " + i);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");

    d = zeros(size(inlierIndicies, 1), 1);

    norm_test = dot(model.Normal', prev_plane_norm);
    for j = 1:size(inlierIndicies, 1)
        d(j) = norm(xyz_P(:, inlierIndicies(j))', plane.Location(j)); 
    end
    idx = find(d > 0.2);
    outlierIndicies = cat(1, outlierIndicies, idx);

    if(abs(norm_test) > (0.09)) 
        break;
    end

    rest_pt_lidar = select(pt_lidar, outlierIndicies);

    pt_lidar = rest_pt_lidar;
    prev_plane_norm = model.Normal;
    ref_vector = circshift(ref_vector, 1); 
    i = i + 1;

end
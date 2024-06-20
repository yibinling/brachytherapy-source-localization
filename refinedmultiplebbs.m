
clearvars
addpath('bin');
ParamSetting;

%% Make measurement - projection
load Datatest.mat % Ground-truth image
NoiseOn = 0;    


%Initial set up
x0 = LAC_380keV;
x0= imresize3(x0,2);
BBs_3D = zeros(size(x0));
DBS = param.DSD - param.DBD; %Distance from Source(center)to BB
detectorX = param.DSD;
source_location = [param.sourcex,param.sourcey,param.sourcez];
shiftx = param.sourcex;
shiftz = param.sourcez;

%Define BB geoometry and location
BBs = [
        struct('cx', (250+DBS)/param.dx , 'cy',(250 -30)/param.dy, 'cz', (86 + 30)/param.dz);             
        struct('cx', (250+DBS)/param.dx , 'cy',(250 +0)/param.dy, 'cz', (86 - 0)/param.dz);   
        struct('cx', (250+DBS)/param.dx , 'cy',(250 -30+40)/param.dy, 'cz', (86 + 30)/param.dz); 
        struct('cx', (250+DBS)/param.dx , 'cy',(250 +30)/param.dy, 'cz', (86- 30)/param.dz);
        struct('cx', (250+DBS)/param.dx , 'cy',(250 +40)/param.dy, 'cz', (86- 0)/param.dz);   
        struct('cx', (250+DBS)/param.dx , 'cy',(250 +30+40)/param.dy, 'cz', (86- 30)/param.dz);
    
];

n = numel(BBs);
bb_centers = zeros(n, 3); 
for i = 1:n
    bb_centers(i, 1) = (BBs(i).cx - 256*2) * param.dx;
    bb_centers(i, 2) = -(BBs(i).cy - 256*2) * param.dy;
    bb_centers(i, 3) = (BBs(i).cz - 86*2) * param.dz;
end
tungsten_attenuation_coefficient = 0.4041 *10;  

for i = 1:numel(BBs)
    [x, y, z] = ndgrid(1:param.nx, 1:param.ny, 1:param.nz);
    distances = sqrt((x - BBs(i).cx).^2 + (y - BBs(i).cy).^2 + (z - BBs(i).cz).^2);
    BBs_3D(distances <=5 /param.dx) = tungsten_attenuation_coefficient;
end

x1 = x0 + BBs_3D;
projWithoutBB = CTprojection2(x0,param);
proj = CTprojection2(x1,param);
maxPixelvalue = max(max(projWithoutBB));

save proj.mat proj 
proj0 = proj(:, :, 1);
figure;
imshow(proj0, []);

%Binarize projection image with threshold = maxPixelvalue of projWithoutBB
bw =  imbinarize(proj0,maxPixelvalue);
minSize = 15;
bw= bwareaopen(bw,minSize);
figure;
imshow(bw)
title("Circular Objects in Binary Value")

stats = regionprops(bw, 'Area', 'PixelIdxList');
maxSize = 6000;

% Remove large objects
for k = 1:length(stats)
    if stats(k).Area > maxSize
        bw(stats(k).PixelIdxList) = 0;  
    end
end

% Fill any gaps in circular objects
se = strel("disk",5);
bw = imclose(bw,se);
imshow(bw)

% Display the label matrix and draw each boundary.
[B,L] = bwboundaries(bw,"noholes");
figure;
imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2),boundary(:,1),"r",LineWidth=2)
end

% Obtain circularity of each boundary. If the circularity exceeds the threshold, calculate the position of the centroid 
centroids = [];
stats2 = regionprops(L,"Circularity","Centroid");
threshold = 0.7;
for k = 1:length(B)
  boundary = B{k};
  % sprintf("%2.2f",circ_value);
  if stats2(k).Circularity > threshold
    centroid = stats2(k).Centroid;
    centroids = [centroids; stats2(k).Centroid];
    plot(centroid(1),centroid(2),"ko");
  end

  text(boundary(1,2)-35,boundary(1,1)+13,sprintf("%2.2f",stats2(k).Circularity),Color="y",...
       FontSize=14,FontWeight="bold")

end
title("Centroids of Circular Objects and Circularity Values")


% Sort the centers from left to right
centers = sortrows(centroids, 1);

%Plot the projected projected BB location
p = zeros(n, 3); 
for i = 1:n
    p(i,:) = [param.DSD,250-centers(i,1)*param.du,250-centers(i,2)*param.dv];
end

% Define the 3D grid for the voxel space
[x, y, z] = meshgrid(param.xs, param.ys, param.zs);
figure;
slice(x, y, z, BBs_3D, size(x,2)/2, size(y,2)/2, size(z,2)/2);
hold on;
xlabel('x-axis (mm)');
ylabel('y-axis (mm)');
zlabel('z-axis (mm)');
title('3D Voxel Space Visualization');
shading interp;

%Visualize the detector
detectorLength = 500; 
detectorWidth = 500; 
[yPlane, zPlane] = meshgrid(linspace(-detectorLength/2, detectorLength/2, 2), linspace(-detectorWidth/2, detectorWidth/2, 2));
xPlane = repmat(detectorX, size(yPlane));
mesh(xPlane, yPlane, zPlane, 'FaceColor', 'none', 'EdgeColor', 'black', 'LineStyle', '--');
hold on;

%Plot the BB locations, source location, and projected BB locations
scatter3(bb_centers(:, 1), bb_centers(:, 2), bb_centers(:, 3),'k', 'filled'); %BBs

scatter3(source_location(1), source_location(2), source_location(3),100,'red','filled') %Source
trueSourceLabel = sprintf('Source: (%.2f, %.2f, %.2f)', source_location(1), source_location(2), source_location(3));
text(source_location(1), source_location(2), source_location(3), trueSourceLabel, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

hold on; 
colors = ['m','b','y','g','c','r'];
for i = 1:n  %Projected BBs
    scatter3(p(i,1), p(i,2), p(i,3), colors(i),'filled');
end

%Plot the extensions from projected source location to the BB locations
extension_length = 300;
for i = 1:n %3D
    dx = p(i,1) - bb_centers(i,1);
    dy = p(i,2) - bb_centers(i,2);
    dz = p(i,3) - bb_centers(i,3);

    len = sqrt(dx^2 + dy^2 + dz^2);
    dx = dx / len;
    dy = dy / len;
    dz = dz / len;

    extended_xi = [bb_centers(i,1) - extension_length*dx, bb_centers(i,1) + extension_length*dx];
    extended_yi = [bb_centers(i,2) - extension_length*dy, bb_centers(i,2) + extension_length*dy];
    extended_zi = [bb_centers(i,3) - extension_length*dz, bb_centers(i,3) + extension_length*dz];

    plot3(extended_xi, extended_yi, extended_zi, colors(i), 'LineStyle', '--');
    hold on;
end

%Plot the projected source location

projected_source_location = findOptimalPoint(bb_centers, p,param);

scatter3(projected_source_location(1),projected_source_location(2),projected_source_location(3), 100, 'k','filled')
projectedSourceLabel = sprintf('Projected Source: (%.2f, %.2f, %.2f)', projected_source_location(1), projected_source_location(2), projected_source_location(3));
text(projected_source_location(1), projected_source_location(2), projected_source_location(3), projectedSourceLabel, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
hold off;

%Output the deviation
deviation_source = norm(projected_source_location - source_location);
disp('Deviation from True Source location:');
disp(deviation_source);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculating line-detector intersection%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intersections = zeros(n, 3); 
deviations = zeros(n, 1); 

for i = 1:n
    bbX = bb_centers(i,1);
    bbY = bb_centers(i,2);
    bbZ = bb_centers(i,3);

    dirVec = [bbX - source_location(1), bbY - source_location(2), bbZ - source_location(3)];
    t = (detectorX - source_location(1)) / dirVec(1);

    intersectY = source_location(2) + t * dirVec(2);
    intersectZ = source_location(3) + t * dirVec(3);

    intersections(i,:) = [detectorX, intersectY, intersectZ]; % Store intersections

    % Calculate deviations between calculated intersections and projected locations in 'p'
    projY = p(i,2);
    projZ = p(i,3);

    deviation = sqrt((projY - intersectY)^2 + (projZ - intersectZ)^2);
    deviations(i) = deviation;
end

% Output the deviations
disp('Deviations between calculated intersections and projected BB locations:');
disp(deviations);

% % Visualization
figure;
[x, y, z] = meshgrid(param.xs, param.ys, param.zs);
slice(x, y, z, BBs_3D, size(x,2)/2, size(y,2)/2, size(z,2)/2);
hold on;
xlabel('x-axis (mm)');
ylabel('y-axis (mm)');
zlabel('z-axis (mm)');
title('3D Visualization of source-to-BB line intersection on the Detector');
shading interp;
%detector plane
detectorLength = 500; 
detectorWidth = 500;
[yPlane, zPlane] = meshgrid(linspace(-detectorLength/2, detectorLength/2, 2), linspace(-detectorWidth/2, detectorWidth/2, 2));
xPlane = repmat(detectorX, size(yPlane));
mesh(xPlane, yPlane, zPlane, 'FaceColor', 'none', 'EdgeColor', 'black', 'LineStyle', '--');

%Plot the BB locations, source location, and projected BB locations
scatter3(source_location(1), source_location(2), source_location(3), 100, 'red', 'filled');
text(source_location(1), source_location(2), source_location(3), ' Source', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

scatter3(bb_centers(:,1), bb_centers(:,2), bb_centers(:,3), 'k', 'filled');

for i = 1:n  %Projected BBs
    scatter3(p(i,1), p(i,2), p(i,3), colors(i),'filled');
end
hold on;

% Draw lines from the source through BBs to their intersection points on the detector
for i = 1:n
    line([source_location(1), intersections(i,1)], [source_location(2), intersections(i,2)], [source_location(3), intersections(i,3)], 'Color', colors(i), 'LineStyle', '--');
end

% Plot the calculated intersection points on the detector plane
scatter3(intersections(:,1), intersections(:,2), intersections(:,3), 'k', 'filled');

grid on;
axis equal;
hold off;

function optimal_point = findOptimalPoint(bb_centers, p, param) %3D
    n = size(bb_centers, 1); 

    function d = totalDistanceSquared(x)
        d = 0;
        for i = 1:n
            linePoint = bb_centers(i,:);
            lineDir = p(i,:) - bb_centers(i,:);
            lineDir = lineDir / norm(lineDir); 

            vec = x - linePoint;
            distVec = vec - dot(vec, lineDir) * lineDir;
            d = d + norm(distVec)^2;
        end
    end 

    x0 = [param.sourcex,param.sourcey,param.sourcez]; % Initial guess
    optimal_point = fminunc(@totalDistanceSquared, x0);
end


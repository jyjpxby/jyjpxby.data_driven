clearvars
close all
%% Define parameter
data1 = readtable('reaction_Setting.xlsx');
plie = table2array(data1(1,1)); % Number of classes
sheetnum = table2array(data1(1,2)); % Read the Excel sheet number
dataen = data1{1,4}; % The name of the file from which the imported data is extracted
delimiter = ' ';  % Set separator
excelname = strjoin(dataen, delimiter);
[~, sheetNames] = xlsfinfo(excelname);
for i = 1:numel(sheetNames)
    sheetData = xlsread(excelname, sheetNames{i});
    ss(1,i) = {sheetData};
end % The data for each worksheet is read through a loop
A = cell2mat(ss(1,sheetnum));
r = 10; % Base radius
h = 5; % Height/Number of layers
row = table2array(data1(1,3)); %Extract vertex
pp = 1:1:plie;
%% data importing
for i = 1:size(A,1) 
    Ad(i,1) = A(i,2); % Extract data field
end
peak = Ad(row, :);
n = plie;
nn = 8*nchoosek(n,2); % Number of points per layer
combinations = nchoosek(1:n, 2);
combinations_column = reshape(combinations.', [], 1);% Data fixed point
combinations_column = combinations_column';
combinations_column11 = combinations_column;
for i = 1:size(combinations_column11,2)
    if combinations_column11(1,i)>=row
        combinations_column11(1,i)=combinations_column11(1,i)+1;
    end
end
newAd = Ad;
newAd(row, :) = [];
num_layers = h-1; % number of plies
plies = nn;% The number of points per layer
increment = h / (num_layers+1); % Calculate the height increment for each layer
layer_heights = increment:increment:h; % The height of each floor
x_points = []; % Stores the X-coordinates of points on each layer
y_points = []; % Stores the Y-coordinates of points on each layer
z_points = []; % Stores the Z-coordinates of points on each layer
for i = 1:num_layers
    radius = (h - layer_heights(i)) / h * r; % The radius of the corresponding layer
    num_points = nn; % The number of points per layer
    theta = linspace(0, 2*pi, num_points + 1); % polar angle
    theta = theta(1:end-1);
    x = radius * cos(theta); % x-coordinate
    y = radius * sin(theta); % y-coordinate
    z = layer_heights(i) * ones(size(theta)); % z-coordinate
    x_points = [x_points, x];
    y_points = [y_points, y];
    z_points = [z_points, z];
end
%% coning
for i = 1:num_layers
    x = [0, x_points]; % Add vertex
    y = [0, y_points];
    z = [h, z_points];
    %plot3(x, y, z, 'ro');
end
nlie = size(x_points, 2);
A2 = zeros(1, nlie);
A2(1,1) = peak;
for i = 1:nn/4
    for ii = 1:plie
        if combinations_column(1,i) == pp(1,ii)
        A2(1,2+4*(i-1)) = newAd(ii,1);
        jyj(1,i) = newAd(ii,1);
        end
    end
end % Add endpoint value
for i = 1:nn/4
    if i ~= plie
        A2(1,4+4*(i-1)) = (A2(1,2+4*(i-1))+A2(1,2+4*(i)))/2;
    else
        A2(1,4+4*(i-1)) = (A2(1,2+4*(i-1))+A2(1,2))/2;
    end
end % Adds endpoint midpoint values
endmidpoint = [];
endmidpointx = [];
endmidpointy = [];
for i = 1:nn/4
    endmidpoint = [endmidpoint,A2(1,4+4*(i-1))]; % Record the endpoint midpoint
    endmidpointx = [endmidpointx,x(1,4+4*(i-1))];
    endmidpointy = [endmidpointy,y(1,4+4*(i-1))];
end
[~,~,idz] = unique(endmidpoint);
idzup = ones(size(idz,1),1)*size(idz,1);
idzdown = zeros(size(idz,1),1);
for i = 1:size(idz,1)
    if endmidpoint(1,i) <= peak
        idzdown(i,1) = idz(i,1);
    end
end
for i = 1:size(idz,1)
    if endmidpoint(1,i) > peak
        idzup(i,1) = idz(i,1);
    end
end
for i = 1:nn/2
    if i ~= plie*2
        A2(1,3+2*(i-1)) = (A2(1,2+2*(i-1))+A2(1,2+2*(i)))/2;
    else
        A2(1,3+2*(i-1)) = (A2(1,2+2*(i-1))+A2(1,2))/2;
    end
end % Adds the midpoint value for the midpoint of the endpoint
for i = 1:num_points
        A2(1,i+1+2*num_points) = (A2(1,i+1)+A2(1,1))/2;
end % Adds a boundary midpoint value
for i = 1:num_points
        A2(1,i+1+3*num_points) = (A2(1,i+1+2*num_points)+A2(1,1))/2;
end % Add a boundary by the inside 1/4 value
for i = 1:num_points
        A2(1,i+1+num_points) = (A2(1,i+1+2*num_points)+A2(1,i+1))/2;
end % Add a boundary near the outer 3/4 value
theta = linspace(0, 2*pi, nn);
X = r * cos(theta);
Y = r * sin(theta);
radius = r;
thetaa = linspace(0, 2*pi, num_points + 1); % polar angle
thetaa = thetaa(1:end-1);
xout = radius * cos(thetaa); % x-coordinate
yout = radius * sin(thetaa); % y-coordinate
zout = ones(1,nn); % z-coordinate
zout = (-1)*zout;
x_points = [x_points, xout];
y_points = [y_points, yout];
z_points = [z_points, zout];
x_bottom = [0,x_points];
y_bottom = [0,y_points];
z_bottom = [h,z_points];
zz_bottom = [A2,zout];
%% interpolate
interp_method = 'linear'; % interpolation method
F1 = scatteredInterpolant(x', y', z');[xq1, yq1] = meshgrid(-r:0.2:r, -r:0.2:r);
[xq1, yq1] = meshgrid(linspace(min(X(:)), max(X(:)), 200), linspace(min(Y(:)), max(Y(:)), 200));
zq1 = griddata(x_bottom(:), y_bottom(:), z_bottom(:), xq1, yq1, interp_method); % Conical interpolation
zq2 = griddata(x_bottom(:), y_bottom(:), zz_bottom(:), xq1, yq1, interp_method);
radius = r;
theta2 = linspace(0, 2*pi, (nn/4) + 1); % polar angle
theta2 = theta2(1:end-1);
xout_point = radius * cos(theta2); % x-coordinate
yout_point = radius * sin(theta2); % y-coordinate
zout_point = ones(1,nn/4); % z-coordinate
zout_point = (-1)*zout_point; % Data point location
%% Coordinate display position
radius2 = r+1.5;
theta3 = linspace(0, 2*pi, (nn/4) + 1); % polar angle
theta3 = theta3(1:end-1);
xout_point2 = radius2 * cos(theta3); % x-coordinate
yout_point2 = radius2 * sin(theta3); % y-coordinate
%% Draw surface
%figure;surf(xq1, yq1, zq1,'FaceAlpha',0.8);xlabel('x');ylabel('y');zlabel('z');title('The interpolation surface of the cone');
min_bl = min(Ad);
min_bl = abs(min_bl);
figure;
hold on
s1 = surf(xq1, yq1, zq2,'FaceAlpha',0.6);
xlabel('x');
ylabel('y');
zlabel('z');
title('Interpolating surface of the plane');
azimuth1 = 0; % orientation
elevation1 = 90; % pitch angle
rotation1 = 0; % rotation angle
view(azimuth1, elevation1)
camroll(rotation1)
s1.EdgeColor = 'none';
markerSize = 10;  % Point size
markerColor = 'k';  % Point color
markerStyle = 'p';  % Point shape
plot3(xout_point, yout_point, zout_point, [markerStyle, markerColor], 'MarkerSize', markerSize)
for i4 = 1:numel(xout_point2)
    text1 = ['(' num2str(combinations_column11(1,i4)) ')'];
    text2 = {num2str(idz(i4,1))};
    text(xout_point2(i4), yout_point2(i4), text1, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    peakxyz = [0,0,peak];
    endpoint = [xout_point(i4), yout_point(i4),0];
    line1 = [peakxyz;endpoint];
    if jyj(1,i4) > peak
        p1 = plot3(line1(:,1),line1(:,2),line1(:,3),'--','Color','k');
    end
    if jyj(1,i4) <= peak
        p1 = plot3(line1(:,1),line1(:,2),line1(:,3),'--','Color','r');
    end
    p1.LineWidth = 1.5;
    if idz(i4,1) == min(idzup)
        endmid = [endmidpointx(i4), endmidpointy(i4),endmidpoint(1,i4)];
        line2 = [endmid;peakxyz];
        p2 = plot3(line2(:,1),line2(:,2),line2(:,3),'-','Color','k');
        p2.LineWidth = 2;
    end
    if idz(i4,1) == max(idzdown)
        endmid = [endmidpointx(i4), endmidpointy(i4),endmidpoint(1,i4)];
        line2 = [endmid;peakxyz];
        p2 = plot3(line2(:,1),line2(:,2),line2(:,3),'-','Color','r');
        p2.LineWidth = 2;
    end
    %if endmidpoint(1,i4) <= peak
        %text(endmidpointx(i4), endmidpointy(i4), text2,'Color','red', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    %end
    %if endmidpoint(1,i4) > peak
        %text(endmidpointx(i4), endmidpointy(i4), text2,'Color','black', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    %end
end % Add dividers and minimum energy intervals
xlim([-r-3, r+3]);
ylim([-r-3, r+3]);
hold off
zq3 = zq1+zq2/min_bl;
figure;
hold on
s2 = surf(xq1, yq1, zq3,'FaceAlpha',0.8);
xlabel('x');
ylabel('y');
zlabel('z');
title('Interpolation surface of the sum of the squares');
azimuth2 = 45; % orientation
elevation2 = 45; % pitch angle
rotation2 = 0; % rotation angle
view(azimuth2, elevation2)
camroll(rotation2)
s2.EdgeColor = 'none';
markerSize = 10;  % Point size
markerColor = 'k';  % Point color
markerStyle = 'p';  % Point shape
plot3(xout_point, yout_point, zout_point, [markerStyle, markerColor], 'MarkerSize', markerSize)
for i4 = 1:numel(xout_point)
    text1 = ['(' num2str(combinations_column11(1,i4)) ')'];
    text(xout_point2(i4), yout_point2(i4), text1, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    peakxyz = [0,0,h+peak/min_bl];
    endpoint = [xout_point(i4), yout_point(i4),-1];
    line1 = [peakxyz;endpoint];
    p1 = plot3(line1(:,1),line1(:,2),line1(:,3),'--','Color','k');
    p1.LineWidth = 1.5;
    if endmidpoint(1,i4) == min(endmidpoint)
        endmid = [endmidpointx(i4), endmidpointy(i4),endmidpoint(1,i4)/min_bl];
        line2 = [endmid;peakxyz];
        p2 = plot3(line2(:,1),line2(:,2),line2(:,3),'-','Color','b');
        p2.LineWidth = 2;
    end
end % Add dividers and minimum energy intervals
xlim([-r-3, r+3]);
ylim([-r-3, r+3]);
hold off

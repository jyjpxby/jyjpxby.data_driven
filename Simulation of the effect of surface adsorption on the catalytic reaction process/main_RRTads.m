clearvars
close all
tic
%% Define parameter
% Reading Excel files
data1 = readtable('ads_Settings.xlsx');
file = 'ads_Settings.xlsx';
sheet = 'Sheet1';
startColumn = 'G';
endColumn = 'I';
columnRange = [startColumn ':' endColumn];
point_p = xlsread(file, sheet, columnRange); % Mixing ratio of each substance
%% Extract literal data
dataen = data1{1,10}; % The name of the file from which the imported data is extracted
delimiter = ' ';  % Set separator
excelname = strjoin(dataen, delimiter);
plie = table2array(data1(1,1)); % Class quantity
sheetnum = table2array(data1(1,2)); % Read the Excel sheet number
sub = table2array(data1(1,3)); % Vertex number
alpha = 1-table2array(data1(1,4)); % confidence region
r1 = 10; % Base radius

%% Data import reads all worksheets in an Excel file
[~, sheetNames] = xlsfinfo(excelname);
for i = 1:numel(sheetNames)
    sheetData = xlsread(excelname, sheetNames{i});
    ss(1,i) = {sheetData};
end % The data for each worksheet is read through a loop
sd = cell2mat(ss(1,sheetnum));
sdsize = size(sd);
wz = sd(1:3,2:sdsize(1,2));
wzpr = sd(4:6,2:sdsize(1,2));
absd = sd(8:sdsize(1,1),2:sdsize(1,2));
for i = 1:size(absd,2)
    datadd = absd(ceil(size(absd,1)/3*2):size(absd,1),i);
    nv = length(datadd); % Number of data samples
    mean_value = mean(datadd); % Data mean
    std_value = std(datadd); % Data standard deviation
    t_score = tinv(1-alpha/2, nv-1); % T-distributed loci, used to calculate confidence intervals
    margin_error = t_score * std_value / sqrt(nv); % Marginal error
    lower_bound = mean_value - margin_error; % Lower bound of confidence interval
    upper_bound = mean_value + margin_error; % Upper bound of confidence interval
    absda(1,i) = mean_value;
    absdaup(1,i) = upper_bound;
    absdadown(1,i) = lower_bound;
end
pro = [];
for i = 1:size(wz,2)
    for j = 1:size(wz,1)
        if wz(j,i) == sub || wz(j,i) == 0
            for ii = 1:3
                pro(ii,i) = wz(ii,i);
                pro(ii+3,i) = wzpr(ii,i);
                pro(7,i) = absda(1,i);
                pro(8,i) = absdaup(1,i);
                pro(9,i) = absdadown(1,i);
            end
        end
    end
end
column_sum = sum(pro); % Calculate the sum of each column to determine which columns are all zeros
nonzero_columns = (column_sum ~= 0); % Creates a logical vector that sets the element corresponding to an all-zero column to false
pro = pro(:, nonzero_columns); % Use logical indexing and slicing operations to remove all zero columns
for i = 1:size(pro,2)
    for j = 1:size(wz,1)
        if pro(j,i) == 0
            temp1 = pro(1, i);
            pro(1, i) = pro(j, i);
            pro(j, i) = temp1;
            temp2 = pro(1+size(wz,1), i);
            pro(1+size(wz,1), i) = pro(j+size(wz,1), i);
            pro(j+size(wz,1), i) = temp2;
        end
    end
end % Substrate sort transposition
for i = 1:size(pro,2)
    for j = 1:size(wz,1)
        if pro(j,i) == sub
            temp1 = pro(1, i);
            pro(1, i) = pro(j, i);
            pro(j, i) = temp1;
            temp2 = pro(1+size(wz,1), i);
            pro(1+size(wz,1), i) = pro(j+size(wz,1), i);
            pro(j+size(wz,1), i) = temp2;
        end
    end
end % Substrate sort transposition
for i = 1:size(pro,2)
    for j = 2:size(wz,1)
        if pro(j,i) > sub
            pro(j,i) = pro(j,i)-1;
        end
    end
end % Product number
Ad = pro(7,:);
n_n = size(point_p,1);
jyj1 = point_p(:,1);
jyj1 = jyj1';
jyj2 = point_p(:,2);
jyj2 = jyj2';
jyj3 = point_p(:,3);
jyj3 = jyj3';
point_n = size(point_p,1);
hh = [];
hg = [];
for j = 1:1:point_n % number of plies
    h_n = 1-(jyj1(1,j)/(jyj1(1,j)+jyj2(1,j)+jyj3(1,j)));
    hh = [hh,h_n];
end
hh1 = sort(hh);
h3 = unique(hh1);
h3_counts = histc(hh1, h3);
h = size(h3,2);% Height/Number of layers
for jj = 1:1:point_n % angular deviation
    h_m = jyj3(1,jj)/(jyj2(1,jj)+jyj3(1,jj));
    hg = [hg,h_m];
end
hg1 = sort(hg);
h4 = unique(hg1);
h4_counts = histc(hg1, h4);
gh = size(h4,2);
for i=1:size(pro,2)
    if pro(5,i)==0 && pro(6,i)==0
        peak = pro(7,i);
        peakup = pro(8,i);
        peakdown = pro(9,i);
    end
end%Extract vertex
n = plie;
nm = 2*nchoosek(n,2); % The number of outermost points
combinations = nchoosek(1:n, 2);
combinations_column = reshape(combinations.', [], 1);% Data fixed point
combinations_column = combinations_column';
combinations_column11 = combinations_column;
for i = 1:size(combinations_column11,2)
    if combinations_column11(1,i)>=sub
        combinations_column11(1,i)=combinations_column11(1,i)+1;
    end
end
%% Calculate the points on each layer of the cone
layer_heights = [];
for i = 1:n_n
    layer_heights(1,i) = h*hh(1,i);
end
theta = linspace(0, 2*pi, nm+1); % polar angle
theta = theta(1:end-1); % Calculate the starting Angle of the points in each group
dtheta = theta(1,2)-theta(1,1); % Calculate the Angle difference between each set of points
x3 = [];y3 = [];z3 = [];z5 = [];rd = [];td = [];
lh = sort(layer_heights);
h5 = unique(lh);
h5_counts = histc(lh, h5);
lhn = size(h5,2);
for i = 1:nm
    for ii = 1:n_n
        radius = hh(ii) * r1; % The radius of the corresponding point
        thetan = (hg(1,ii)*dtheta)+theta(1,i); % The deflection Angle of the corresponding point
        x = radius * cos(thetan); % x-coordinate
        y = radius * sin(thetan); % y-coordinate
        rd(ii,i) = radius;
        td(ii,i) = thetan;
        z = (r1-radius)/r1*h; % z-coordinate
        x3 = [x3,x];
        y3 = [y3,y];
        z3 = [z3,z];
    end
end
x = [0, x3]; % Add conical vertex
y = [0, y3];
z = [h, z3]; % Conical point
Aup3(1,1) = peakup; % Add data vertex
Adown3(1,1) = peakdown; % Add data vertex
A3(1,1) = peak; % Add data vertex
z4 = [];zup = [];zdown = [];
for i = 1:nm
    if i <= nm-1
        for j = 1:size(pro,2)
            if ((combinations_column(1,i)==pro(2,j)) && (combinations_column(1,i+1)==pro(3,j))) || ((combinations_column(1,i)==pro(2,j)) && (pro(3,j)==0))
                for i1 = 1:nm
                    for i2 = 1:n_n
                        if ((((hg(1,i2)*dtheta)+theta(1,i1)) == ((pro(5,j)/(pro(5,j)+pro(6,j))*dtheta)+theta(1,i)))) && (hh(1,i2) == (1-(pro(4,j)/(pro(4,j)+pro(5,j)+pro(6,j)))))
                            z4(i1,i2) = pro(7,j);
                            zup(i1,i2) = pro(8,j);
                            zdown(i1,i2) = pro(9,j);
                        end
                    end
                end
            end
            if ((combinations_column(1,i)==pro(3,j)) && (combinations_column(1,i+1)==pro(2,j))) || ((combinations_column(1,i)==pro(2,j)) && (pro(3,j)==0))
                for i1 = 1:nm
                    for i2 = 1:n_n
                        if ((((hg(1,i2)*dtheta)+theta(1,i1)) == ((pro(6,j)/(pro(5,j)+pro(6,j))*dtheta)+theta(1,i)))) && (hh(1,i2) == (1-(pro(4,j)/(pro(4,j)+pro(5,j)+pro(6,j)))))
                            z4(i1,i2) = pro(7,j);
                            zup(i1,i2) = pro(8,j);
                            zdown(i1,i2) = pro(9,j);
                        end
                    end
                end
            end
        end
    end
    if i == nm
        for j = 1:size(pro,2)
            if ((combinations_column(1,i)==pro(2,j)) && (combinations_column(1,1)==pro(3,j))) || ((combinations_column(1,i)==pro(2,j)) && (pro(3,j)==0))
                for i1 = 1:nm
                    for i2 = 1:n_n
                        if (((hg(1,i2)*dtheta)+theta(1,i1) == (pro(5,j)/(pro(5,j)+pro(6,j))*dtheta)+theta(1,i))) && (hh(1,i2) == 1-(pro(4,j)/(pro(4,j)+pro(5,j)+pro(6,j))))
                            z4(i1,i2) = pro(7,j);
                            zup(i1,i2) = pro(8,j);
                            zdown(i1,i2) = pro(9,j);
                        end
                    end
                end
            end
            if ((combinations_column(1,i)==pro(3,j)) && (combinations_column(1,1)==pro(2,j))) || ((combinations_column(1,i)==pro(2,j)) && (pro(3,j)==0))
                for i1 = 1:nm
                    for i2 = 1:n_n
                        if (((hg(1,i2)*dtheta)+theta(1,i1) == (pro(6,j)/(pro(5,j)+pro(6,j))*dtheta)+theta(1,i))) && (hh(1,i2) == 1-(pro(4,j)/(pro(4,j)+pro(5,j)+pro(6,j))))
                            z4(i1,i2) = pro(7,j);
                            zup(i1,i2) = pro(8,j);
                            zdown(i1,i2) = pro(9,j);
                        end
                    end
                end
            end
        end
    end
end% Add endpoint values and data point values
radius = r1+1;
radius1 = r1+2.5;
thetaa = linspace(0, 2*pi, nm + 1); % polar angle
thetaa = thetaa(1:end-1);
xout = radius * cos(thetaa); % x-coordinate
yout = radius * sin(thetaa); % y-coordinate
xout1 = radius1 * cos(thetaa); % x-coordinate
yout1 = radius1 * sin(thetaa); % y-coordinate
zout0 = ones(1,nm); % z-coordinate
zoutb = mean(pro(7,:));
zout2 = (zoutb)*zout0;
zout = (-1)*zout0;zoutup = mean(pro(8,:))*zout0;zoutdown = mean(pro(9,:))*zout0;
x3 = [x3, xout];
y3 = [y3, yout];
z3 = [z3, zout];
z5 = [z5,zout2];
z4 = reshape(z4,1,[]);zup = reshape(zup,1,[]);zdown = reshape(zdown,1,[]);
A3 = [A3,z4];Aup3 = [Aup3,zup];Adown3 = [Adown3,zdown];
A3 = [A3,zout2];Aup3 = [Aup3,zoutup];Adown3 = [Adown3,zoutdown];
x_bottom = [0,x3];
y_bottom = [0,y3];
z_bottom = [h,z3];
zz_bottom = A3;zup_bottom = Aup3;zdown_bottom = Adown3;
%% interpolate
interp_method = 'linear'; % interpolation method
[xq1, yq1] = meshgrid(linspace(min(x_bottom(:)), max(x_bottom(:)), 200), linspace(min(y_bottom(:)), max(y_bottom(:)), 200));
zq1 = griddata(x_bottom(:), y_bottom(:), z_bottom(:), xq1, yq1, interp_method); % Conical interpolation
zqup1 = griddata(x_bottom(:), y_bottom(:), zup_bottom(:), xq1, yq1, interp_method); % Upper interpolation
zqdown1 = griddata(x_bottom(:), y_bottom(:), zdown_bottom(:), xq1, yq1, interp_method); % Lower interpolation
zq2 = griddata(x_bottom(:), y_bottom(:), zz_bottom(:), xq1, yq1, interp_method); % Plane interpolation
min_bl = min(Ad);
min_bl = abs(min_bl);
zq3 = zq1+zq2/min_bl;
for i = 1:200
    for j = 1:200
        dzq(i,j) = abs(zqup1(i,j)-zqdown1(i,j));
    end
end
ddzq = max(dzq);
zqup1 = 0.3+2*ddzq+zqup1;
zqdown1 = -0.3-2*ddzq+zqdown1;
%% Drawing
figure(1)
axis([-5-r1 5+r1 -5-r1 5+r1 min(Adown3)-5 max(Aup3)+5])
surf(xq1,yq1,zqup1,'FaceAlpha',0.2,EdgeColor = 'none')
hold on
azimuth3 = 0; % orientation
elevation3 = 90; % pitch angle
rotation3 = 0; % rotation angle
view(azimuth3, elevation3)
camroll(rotation3)
markerSize = 10;  % Point size
markerColor = 'k';  % Point color
markerStyle = 'p';  % Point shape
plot3(xout, yout, zout2, [markerStyle, markerColor], 'MarkerSize', markerSize)
for i4 = 1:numel(xout1)
    text1 = ['(' num2str(combinations_column11(1,i4)) ')'];
    text(xout1(i4), yout1(i4),zout2(i4), text1 , 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    peakxyz = [0,0,peak];
    endpoint = [xout(i4), yout(i4),zoutb];
    line1 = [peakxyz;endpoint];
    p1 = plot3(line1(:,1),line1(:,2),line1(:,3),'--','Color','k');
    p1.LineWidth = 1.5;
end % Add a divider
hold on
surf(xq1,yq1,zqdown1,'FaceAlpha',1,EdgeColor = 'none') %Remove the leading % to show the upper limit of energy
hold on
X1 = xq1;Y1 = yq1;Z1 = zqdown1;
XA=X1;YA=Y1;ZA=Z1;
xa=XA';XA=xa(:);
ya=YA';YA=ya(:);
za=ZA';ZA=za(:);
X2 = xq1;Y2 = yq1;Z2 = zqup1;
XB=X2;YB=Y2;ZB=Z2;
xb=XB';XB=xb(:);
yb=YB';YB=yb(:);
zb=ZB';ZB=zb(:);
XX=[XA;XB];
YY=[YA;YB];
ZZ=[ZA;ZB];
obstacle = [XX,YY,ZZ];
x_min = -radius;
y_min = -radius;
z_min = min(pro(9,:));
x_max = radius;
y_max = radius;
z_max = max(pro(8,:));
EPS = table2array(data1(1,5));
numNodes = table2array(data1(1,6));        
q_start.coord = [0 0 peak];
q_start.cost = 0;
q_start.parent = 0.5;
q_goal.cost = 0;
nodes(1) = q_start;
%% Search for possible response paths
for i = 1:1:numNodes
    axqr = rand(1);
    ayqr = rand(1);
    bxqr = rand(1);
    byqr = rand(1);
    xqr = axqr*(r1+2)*cos(bxqr*2*pi);
    yqr = ayqr*(r1+2)*sin(byqr*2*pi);
    zqr = z_min+rand(1)*(z_max-z_min);
    q_rand = [xqr yqr zqr];
    xyii = rem(i+nm-1,nm);
    xyii = xyii+1;
    q_goal.coord = [xout(1,xyii) yout(1,xyii) peak];
    %plot3(q_rand(1), q_rand(2), q_rand(3), 'x', 'Color',  [0 0.4470 0.7410]) % Display random point
    % Break if goal node is already reached
    for j = 1:1:length(nodes)
        if nodes(j).coord == q_goal.coord
            break
        end
    end
    % Pick the closest node from existing list to branch out from
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist_3d(n.coord, q_rand);
        ndist = [ndist tmp];
    end
    [val, idx] = min(ndist);
    q_near = nodes(idx);
    rr = ones(402,402);
    rr = rr*0.4;
    q_new.coord = steer3d(q_rand, q_near.coord, val, EPS);
    if checkPath3(q_rand, q_near.coord, obstacle,rr)
        line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], [q_near.coord(3), q_new.coord(3)], 'Color', 'k', 'LineWidth', 3);
        drawnow
        hold on
        q_new.cost = dist_3d(q_new.coord, q_near.coord) + q_near.cost;
        % Within a radius of r, find all existing nodes
        q_nearest = [];
        r = 0.5;
        neighbor_count = 1;
        for j = 1:1:length(nodes)
            if checkPath3(nodes(j).coord, q_new.coord, obstacle,rr) && (dist_3d(nodes(j).coord, q_new.coord)) <= r
                q_nearest(neighbor_count).coord = nodes(j).coord;
                q_nearest(neighbor_count).cost = nodes(j).cost;
                neighbor_count = neighbor_count+1;
            end
        end
        % Initialize cost to currently known value
        q_min = q_near;
        C_min = q_new.cost;
        % Iterate through all nearest neighbors to find alternate lower
        % cost paths
        for k = 1:1:length(q_nearest)
            if checkPath3(q_nearest(k).coord, q_new.coord, obstacle,rr) && q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord) < C_min
                q_min = q_nearest(k);
                C_min = q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord);
                line([q_min.coord(1), q_new.coord(1)], [q_min.coord(2), q_new.coord(2)], [q_min.coord(3), q_new.coord(3)], 'Color', 'g');
                hold on
            end
        end
        % Update parent to least cost-from node
        for j = 1:1:length(nodes)
            if nodes(j).coord == q_min.coord
                q_new.parent = j;
            end
        end
        % Append to nodes
        nodes = [nodes q_new];
    end
end
D = [];
for j = 1:1:length(nodes)
    tmpdist = dist_3d(nodes(j).coord, q_goal.coord);
    D = [D tmpdist];
end
disp('All simulated paths have been completed')
hold off
toc
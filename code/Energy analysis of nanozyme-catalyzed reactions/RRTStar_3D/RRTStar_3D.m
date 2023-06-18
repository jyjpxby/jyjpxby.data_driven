clearvars
close all
%% Read data

data1=importdata('steer3d_Settings.txt');%Each time you make a change, replace the value in this file with the corresponding value in Table S4.
data2=data1.data;
x_min = data2(1,1);
y_min = data2(1,2);
z_min = data2(1,3);
x_max = data2(1,4)-x_min;
y_max = data2(1,5)-y_min;
z_max = data2(1,6)-z_min;
A=load('steer3d_Lower_limit.txt');%Each time you make a change, replace the value in this file with the corresponding value in Table S4.
x1=A(:,1);y1=A(:,2);z1=A(:,3);
x1=x1-x_min;
y1=y1-y_min;
z1=z1-z_min;
[X1,Y1,Z1]=griddata(x1,y1,z1,linspace(min(x1),max(x1),201)',linspace(min(y1),max(y1),201),'v4');%interpolation
XA=X1;YA=Y1;ZA=Z1;
xa=XA';
XA=xa(:);
ya=YA';
YA=ya(:);
za=ZA';
ZA=za(:);
B=load('steer3d_Upper_limit.txt');%Each time you make a change, replace the value in this file with the corresponding value in Table S4.
x2=B(:,1);y2=B(:,2);z2=B(:,3);
x2=x2-x_min;
y2=y2-y_min;
z2=z2-z_min;
[X2,Y2,Z2]=griddata(x2,y2,z2,linspace(min(x2),max(x2),201)',linspace(min(y2),max(y2),201),'v4');%interpolation
XB=X2;YB=Y2;ZB=Z2;
xb=XB';
XB=xb(:);
yb=YB';
YB=yb(:);
zb=ZB';
ZB=zb(:);
XX=[XA;XB];
YY=[YA;YB];
ZZ=[ZA;ZB];
obstacle = [XX,YY,ZZ];
W=[200	198.7622486	500
231.6760636	400	500
200	301.9082001	500
168.3239364	400	500
102.5115325	333.7802991	500
400	230.6343476	500
297.4884675	230.6343476	500
380.4202973	291.2841671	500
297.4884675	333.7802991	500
297.4884675	7.32896E-14	500
357.7652486	147.1892728	500
260.2767811	115.3171738	500
357.7652486	44.0433213	500
42.23475141	44.0433213	500
102.5115325	7.32896E-14	500
200	31.87209902	500
139.7232189	115.3171738	500
19.57970272	291.2841671	500
-7.28388E-14	230.6343476	500
42.23475141	147.1892728	500
102.5115325	230.6343476	500
];%The position of the point corresponding to each substance.
xxx=W(:,1);yyy=W(:,2);zzz=W(:,3);
xxx=xxx-x_min;
yyy=yyy-y_min;
zzz=zzz-z_min;
EPS = data2(1,7);
numNodes = data2(1,8);        
q_start.coord = [data2(1,9)-x_min data2(1,10)-y_min data2(1,11)-z_min];
q_start.cost = 0;
q_start.parent = 0;
q_goal.coord = [data2(1,12)-x_min data2(1,13)-y_min data2(1,14)-z_min];
q_goal.cost = 0;
nodes(1) = q_start;
%% Drawing

figure(1)
axis([0 x_max 0 y_max 0 z_max])
surf(X1,Y1,Z1,EdgeColor = 'none')
hold on
scatter3(xxx,yyy,zzz,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5])%Displays the position of the point corresponding to each substance.
hold on
%surf(X2,Y2,Z2,'FaceAlpha',0.5,EdgeColor = 'none') %Remove the leading % to show the upper limit of energy
hold on
%% Search for possible response paths

for i = 1:1:numNodes
    q_rand = [rand(1)*x_max rand(1)*y_max rand(1)*z_max];
    plot3(q_rand(1), q_rand(2), q_rand(3), 'x', 'Color',  [0 0.4470 0.7410])
    
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
    rr = rr*8;
    q_new.coord = steer3d(q_rand, q_near.coord, val, EPS);
      if checkPath3(q_rand, q_near.coord, obstacle,rr)
          line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], [q_near.coord(3), q_new.coord(3)], 'Color', 'k', 'LineWidth', 4);
          drawnow
          hold on
          q_new.cost = dist_3d(q_new.coord, q_near.coord) + q_near.cost;
    
           % Within a radius of r, find all existing nodes
          q_nearest = [];
          r = 5;
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

%% Search backwards from goal to start to find the optimal least cost path

[val, idx] = min(D);
q_final = nodes(idx);
q_goal.parent = idx;
q_end = q_goal;
nodes = [nodes q_goal];
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)], [q_end.coord(3), nodes(start).coord(3)], 'Color', 'r', 'LineWidth', 4);
    hold on
    q_end = nodes(start);
end


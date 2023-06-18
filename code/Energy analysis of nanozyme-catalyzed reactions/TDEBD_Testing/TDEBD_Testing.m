clearvars
close all
%% Read data

A=load('TDEBD_Testing.txt');
x=A(:,1);y=A(:,2);z=A(:,3);
W=load('TDEBD_P.txt');
xx=W(:,1);yy=W(:,2);zz=W(:,3);
%% Drawing

scatter3(x,y,z)%scatter diagram
figure
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');%interpolation
pcolor(X,Y,Z);
shading interp%False-color graph
figure,contourf(X,Y,Z) %contour map
figure,surf(X,Y,Z);%3D surface
figure,meshc(X,Y,Z)%profile map
view(0,0); 
figure,meshc(X,Y,Z);%3D surface (light color) + contours
hidden off;
figure,contourf(X,Y,Z,20)%Plot the contour lines
[U,V] = gradient(Z,-0.1,-0.1);%Computed gradient
 hold on
 quiver(X,Y,U,V,2,'r')%Gradient diagram
 hold on
 scatter3(xx,yy,zz,'MarkerEdgeColor','k','MarkerFaceColor',[0 .5 .5])
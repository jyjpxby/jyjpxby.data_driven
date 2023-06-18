clearvars
close all
%% Read data

AA = load('Adsorption.txt');
AAA = load('Enthalpy.txt');
BB = AA;
BBB = AAA;
DD =importdata('Settings.txt');
EE =importdata('Removed_paths.txt');
%% Calculate the path matrix

SS = size(EE);
MINN = DD.data(1,2);
MAXX = DD.data(1,3);
MINN1 = DD.data(1,4);
MAXX1 = DD.data(1,5);
BB = BB';
CC = BB - AA;
CC(CC<MINN1) = inf;
CC(CC>MAXX1) = inf;
CC(1:6,1:6) = inf;
CC(7:12,7:12) = inf;
CC(13:18,13:18) = inf;
ii = EE(2,:);
jj = EE(1,:);
CC(sub2ind(size(CC), ii, jj)) = inf;
fid = fopen('Adsorption_path_matrix.txt','w');% Output adsorption and desorption path matrix
fprintf(fid,[repmat('%5.2f\t', 1, size(CC,2)), '\n'], CC');     
fclose(fid);
BBB = BBB';
CCC = BBB - AAA;
CCC(CCC<MINN) = inf;
CCC(CCC>MAXX) = inf;
iii = EE(2,:);
jjj = EE(1,:);
CCC(sub2ind(size(CCC), iii, jjj)) = inf;
fid = fopen('Enthalpy_path_matrix.txt','w');% Output chemical reaction path matrix
fprintf(fid,[repmat('%5.2f\t', 1, size(CCC,2)), '\n'], CCC');     
fclose(fid);
C1 = CCC;
C2 = CC;
C1(C1==inf)=0;
C2(C2==inf)=0;
F1 = C1+C2;
F1(F1==0)=inf;
F1 = abs(F1);
F1(logical(eye(size(F1))))=0;
F1(F1<0) = inf;
F1(F1<MINN) = inf;
F1(F1>MAXX) = inf;
ii = EE(2,:);
jj = EE(1,:);
F1(sub2ind(size(F1), ii, jj)) = inf;
fid = fopen('Path_Matrix.txt','w');% Output total path matrix
fprintf(fid,[repmat('%5.2f\t', 1, size(F1,2)), '\n'], F1');     
fclose(fid);
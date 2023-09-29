clearvars
close all
tic% Timer start
%% Read data
AA =load('Path_Matrix.txt');
data1 = readtable('path_simulation_setting.xlsx');
BB = table2array(data1(1,1));
AAA = table2array(data1(:,11:12));
ssdd = AAA';
ssdd = ssdd(:, ~any(isnan(ssdd)));
names1 = data1(:,6);
names1 = rmmissing(names1);
names1 = cellstr(names1.Active_site_name);
names1 = names1';
%% Simulate feasible paths and search for shortest paths
HH = size(ssdd);
HH = HH(1,2);
Records = [];
i2 = 1;
while i2 <= HH
    A = AA;
    B = 1.5;
    start = ssdd(1,i2);
    dest = ssdd(2,i2);
    s4 = start;
    d4 = dest;
    p = size(A,1);
    Distance = zeros(2,p);
    Distance(1,:) = 1:p;
    Distance(2,1:p) = A(dest,1:p);
    S(1) = dest;
    U = 1:p;
    U(dest) = [];
    new_Distance = Distance;
    D = Distance;
    D(:,dest) = [];
    path = zeros(2,p);
    path(1,:) = 1:p;
    path(2,Distance(2,:)~=inf) = dest;
    while ~isempty(U)
        index = find(D(2,:)==min(D(2,:)),1);
        k = D(1,index);
        S = [S,k];
        U(U==k) = [];
        new_Distance(2,:) = A(k,1:p)+Distance(2,k);
        D = min(Distance,new_Distance);
        path(2,D(2,:)~=Distance(2,:)) = k;
        Distance = D;
        D(:,S) = [];
    end
    dist = Distance(2,start);
    start1 = start;
    start2 = start;
    path1 = path;
    path2 = path;
    while start1 ~= dest
        path3(:,start1) = path2(:,start1);
        next = path1(2,start1);
        if next(1,1)==0
            clearvars -except  AA BB ssdd i2 HH Records names1;
            A =AA;
            B =1.5;
            start = ssdd(2,i2);
            dest = ssdd(1,i2);
            s4 = start;
            d4 = dest;
            p = size(A,1);
            Distance = zeros(2,p);
            Distance(1,:) = 1:p;
            Distance(2,1:p) = A(dest,1:p);
            S(1) = dest;
            U = 1:p;
            U(dest) = [];
            new_Distance = Distance;
            D = Distance;
            D(:,dest) = [];
            path = zeros(2,p);
            path(1,:) = 1:p;
            path(2,Distance(2,:)~=inf) = dest;
            while ~isempty(U)
                index = find(D(2,:)==min(D(2,:)),1);
                k = D(1,index);
                S = [S,k];
                U(U==k) = [];
                new_Distance(2,:) = A(k,1:p)+Distance(2,k);
                D = min(Distance,new_Distance);
                path(2,D(2,:)~=Distance(2,:)) = k;
                Distance = D;
                D(:,S) = [];
            end
            dist = Distance(2,start);
            start1 = start;
            start2 = start;
            path1 = path;
            path2 = path;
            Lia = ismember(inf,Distance);
            if Lia>=1
                clearvars -except  AA BB ssdd i2 HH Records names1;
                path3 = inf;
                close all
                break
            end
        else
            start1 = next;
        end
    end
    if path3~=inf
        path3(:,all(path3==0))=[];
        fprintf('The shortest path found is：');
        path5 = [];
        ii = 1;
        while start2 ~= dest
            fprintf('%d-->',start2);
            next = path(2,start2);
            path5(ii) = start2;
            ii = ii+1;
            start2 = next;
        end
        path5(path5==0)=[];
        path5 = [path5,dest];
        fprintf('%d\n',dest);
        fprintf('The distance corresponding to the shortest path is：%d\n',dist);
        Records(i2,:) = [start,dest,dist];
        save Records.txt -ascii Records
        A(A==inf) = 0;
        [t,s,weights] = find(A);
        names = names1;
        G = digraph(s,t,weights,names); 
        [path4,dist] = shortestpath(G,start,dest);
%% Drawing
        figure;
        hold on
        h = plot(G,'--','EdgeLabel',G.Edges.Weight,'NodeColor','k','NodeFontSize',15,'EdgeColor','#4DBEEE','EdgeFontSize',10);
        h.LineWidth = B;
        highlight(h,path4,'NodeColor','r','EdgeColor','r')
        set(gcf,'Position',[0 0 1456 1232]);
        set(gca,'Position',[.13 .17 .80 .74]);
        set(gca,'Visible','off');
        strl = [num2str(s4),'--',num2str(d4)];
        saveas(h,['Paths\',strl,'.png']);% Change the path of the save folder
        xLimits = xlim;
        yLimits = ylim;
        hold off
        figure;
        hold on
        path6 = (path3-1)*(yLimits(1,2)-yLimits(1,1))/(p-1)+yLimits(1,1);
        xx1 = size(path6);
        xx1 = xx1(1,2);
        xx1 = ones([1,xx1]);
        y1 = path6(1,:);
        y2 = path6(2,:);
        for i=1:BB
            % Set color mapping
            colormap(gca, 'parula');
            % Set the color according to the number of cycles
            colorIndex = mod(i-1, size(colormap, 1)) + 1;% Image Size Adjustment
            fill([xLimits(1,1)-0.03*(xLimits(1,1)-xLimits(1,2)),xLimits(1,1)-0.03*(xLimits(1,1)-xLimits(1,2)),xLimits(1,1),xLimits(1,1)],[(i-1)*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),i*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),i*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),(i-1)*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1)],colorIndex, 'FaceAlpha', 0.6);
            fill([xLimits(1,2),xLimits(1,2),xLimits(1,2)+0.03*(xLimits(1,1)-xLimits(1,2)),xLimits(1,2)+0.03*(xLimits(1,1)-xLimits(1,2))],[(i-1)*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),i*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),i*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1),(i-1)*p/BB*(yLimits(1,2)-yLimits(1,1))/p+yLimits(1,1)],colorIndex, 'FaceAlpha', 0.6);
        end
        plot([xLimits(1,1);xLimits(1,2)],[y1;y2]);
        nnamee = {'start' 'destination'};
        nnamee1 = names1;
        figure_FontSize=26;
        set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
        set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
        set(gca,'xtick',xLimits(1,1):xLimits(1,2)-xLimits(1,1):xLimits(1,2));
        set(gca,'xticklabel',nnamee)
        set(gca,'ytick',yLimits(1,1):(yLimits(1,2)-yLimits(1,1))/(p-1):yLimits(1,2));% Label size adjustment
        set(gca,'yticklabel',nnamee1)
        G1 = digraph(s,t,weights,names1); 
        h1 = plot(G1,'--','EdgeLabel',G1.Edges.Weight,'NodeColor','k','NodeFontSize',15,'EdgeColor','#4DBEEE','EdgeFontSize',10);
        h1.LineWidth = B;
        highlight(h1,path5,'NodeColor','r','EdgeColor','r')
        set(gcf,'Position',[0 0 1456 1232]);
        set(gca,'Position',[.14 .17 .72 .76]);
        set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
        set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
        set(findobj('FontSize',10),'FontSize',figure_FontSize);
        set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3);
        yyaxis right
        set(gca,'ycolor','k');
        set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
        set(gca,'ytick',yLimits(1,1):(yLimits(1,2)-yLimits(1,1))/(p-1):yLimits(1,2));% Label size adjustment
        set(gca,'yticklabel',nnamee1)
        title('unit: eV')
        hold off
        strl = [num2str(s4),'-',num2str(d4)];
        saveas(h1,['Paths and Nodes\',strl,'.png']);% Change the path of the save folder
        i2 = i2+1;
        clearvars -except  AA BB ssdd i2 HH Records names1;
        close all
    else
        i2 = i2+1;
        clearvars -except  AA BB ssdd i2 HH Records names1;
        close all
    end
end
toc% End of time
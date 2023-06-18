clearvars
close all
tic% Timer start
%% Read data

AA =load('Path_Matrix.txt');
BB =importdata('Settings.txt');
ssdd = load('Record_Nodes.txt');
%% Simulate feasible paths and search for shortest paths

HH = size(ssdd);
HH = HH(1,2);
Records = [];
i2 = 1;
while i2 <= HH
    A = AA;
    B = BB;
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
            clearvars -except  AA BB ssdd i2 HH Records;
            A =AA;
            B =BB;
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
                clearvars -except  AA BB ssdd i2 HH Records;
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
        names = {'Co—H2O'  'Co—H2O2'  'Co—O2'  'Co—•OH'  'Co—•OOH'  'Co—O2•-' 'Cu—H2O' 'Cu—H2O2' 'Cu—O2' 'Cu—•OH' 'Cu—•OOH' 'Cu—O2•-' 'Mn—H2O' 'Mn—H2O2' 'Mn—O2' 'Mn—•OH' 'Mn—•OOH' 'Mn—O2•-'};
        G = digraph(s,t,weights,names); 
        [path4,dist] = shortestpath(G,start,dest);
%% Drawing

        h = plot(G,'--','EdgeLabel',G.Edges.Weight,'NodeColor','k','EdgeColor','r');
        h.LineWidth = B.data(1,1);
        highlight(h,path4,'NodeColor','g','EdgeColor','g')
        set(gcf,'Position',[0 0 1456 1232]);
        set(gca,'Position',[.13 .17 .80 .74]);
        set(gca,'Visible','off');
        strl = [num2str(s4),'--',num2str(d4)];
        saveas(h,['Paths\',strl,'.png']);% Change the path of the save folder
        hold off
        figure
        names1 = {'Co—H2O' 'Co—H2O2' 'Co—O2' 'Co—•OH' 'Co—•OOH' 'Co—O2•-' 'Cu—H2O' 'Cu—H2O2' 'Cu—O2' 'Cu—•OH' 'Cu—•OOH' 'Cu—O2•-' 'Mn—H2O' 'Mn—H2O2' 'Mn—O2' 'Mn—•OH' 'Mn—•OOH' 'Mn—O2•-'};
        hold on
        path6 = path3;
        path6 = path6/3;
        path6 = path6-3.3;
        xx1 = size(path6);
        xx1 = xx1(1,2);
        xx1 = ones([1,xx1]);
        x1 = xx1*(-3);
        x2 = xx1*(3.5);
        y1 = path6(1,:);
        y2 = path6(2,:);
        fill([-3,-3,-2.88,-2.88],[-3.15,-1.15,-1.15,-3.15],'g');
        fill([3.38,3.38,3.5,3.5],[-3.15,-1.15,-1.15,-3.15],'g');
        fill([-3,-3,-2.88,-2.88],[-1.15,0.85,0.85,-1.15],'c');
        fill([3.38,3.38,3.5,3.5],[-1.15,0.85,0.85,-1.15],'c');
        fill([-3,-3,-2.88,-2.88],[0.85,2.85,2.85,0.85],'y');
        fill([3.38,3.38,3.5,3.5],[0.85,2.85,2.85,0.85],'y');% Image Size Adjustment
        plot([x1;x2],[y1;y2]);
        nnamee = {'start' 'destination'};
        nnamee1 = {'Co—H2O' 'Co—H2O2' 'Co—O2' 'Co—•OH' 'Co—•OOH' 'Co—O2•-' 'Cu—H2O' 'Cu—H2O2' 'Cu—O2' 'Cu—•OH' 'Cu—•OOH' 'Cu—O2•-' 'Mn—H2O' 'Mn—H2O2' 'Mn—O2' 'Mn—•OH' 'Mn—•OOH' 'Mn—O2•-' '  '};
        figure_FontSize=26;
        set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
        set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
        set(gca,'xtick',-3:6.5:3.5);
        set(gca,'xticklabel',nnamee)
        set(gca,'ytick',-2.97:6/18:3.03);% Label size adjustment
        set(gca,'yticklabel',nnamee1)
        G1 = digraph(s,t,weights,names1); 
        h1 = plot(G1,'--','EdgeLabel',G1.Edges.Weight,'NodeColor','k','EdgeColor','#4DBEEE');
        h1.LineWidth = B.data(1,1);
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
        set(gca,'ytick',-2.15:4.3/17:2.45);% Label size adjustment
        set(gca,'yticklabel',nnamee1)
        title('unit: eV')
        hold off
        strl = [num2str(s4),'-',num2str(d4)];
        saveas(h1,['Paths and Nodes\',strl,'.png']);% Change the path of the save folder
        i2 = i2+1;
        clearvars -except  AA BB ssdd i2 HH Records;
        close all
    else
        i2 = i2+1;
        clearvars -except  AA BB ssdd i2 HH Records;
        close all
    end
end
toc% End of time
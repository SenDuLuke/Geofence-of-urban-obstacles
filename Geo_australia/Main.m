clear
clc
close all
addpath(genpath(pwd));
format long

%% 导入地图
if exist('MelbourneMap.csv','file')~=0
    map=xlsread('MelbourneMap.csv');
else
% 经纬度坐标转换
    data=xlsread('澳大利亚墨尔本局部地区.csv');
    minLatitude=min(data(:,2));
    minLongitude=min(data(:,3));
    [rowm,~]=size(data);
    X=zeros(rowm,2);
    for i=1:rowm
        [X(i,1),X(i,2)]=TransLatLon(data(i,2),data(i,3),0,minLongitude);
    end
%%     
% 转化为栅格数据，5米间隔
    X(:,2)=X(:,2)-min(X(:,2));
    gridm=ceil(max(X(:,2))/6);
    gridn=ceil(max(X(:,1))/6);
    Y=ceil(X/30);
    map=zeros(gridm,gridn);
    
    for i=1:gridm
        for j=1:gridn
            for k=1:rowm
                if (Y(k,1)==j && Y(k,2)==i)
                    map(i,j)=max(map(i,j),data(k,4));
                end
            end
        end
    end
    surf(map);
    % csvwrite('MelbourneMap.csv',map); %地形数据矩阵
end

%% 等高线图
step=10; %等高线步长
contourStep=step:step:ceil(max(max(map))/step)*step;
visible=figure;
contourData=contour(map,contourStep);  %等高线图（contourData表示等高线值及其坐标数量+若干坐标）
%pcolor(X,Y,Z);shading interp %伪彩色图
%set(visible,'visible','off');

%% 等高线处理
if exist('CylinderApproximationNew.csv','file')~=0
    parameters=xlsread('CylinderApproximationNew.csv');
else
    contourBorder=(rem(contourData(1,:),step)==0) & (rem(contourData(2,:),1)==0);   %寻找等高线值（分界点）
    [~,contourBorderPosition]=find(contourBorder==1);                               %分界点在contourData中的位置
    contourData_num=contourData(:,contourBorderPosition);                           %分界点对应的坐标数量
    %向量前后两两元素差值，剔除第一步没能筛选掉的结果
    for i=2:size(contourData_num,2)
        if(contourData_num(1,i)-contourData_num(1,i-1)>step)
            contourBorder(contourBorderPosition(i))=0;
        end
    end
    %重新获取一遍
    [~,contourBorderPosition]=find(contourBorder==1);
    contourData_num=contourData(:,contourBorderPosition);
    
    %%提取单个等高线坐标
    parameters=zeros(size(contourData_num,2),4); %参数汇总（高度、半径、中心坐标）
    for i=1:size(contourData_num,2)
        contourBorderOne=zeros(2,contourData_num(2,i));
        contourBorderOne(1,:)=contourData(1,(contourBorderPosition(i)+1):(contourBorderPosition(i)+contourData_num(2,i)));
        contourBorderOne(2,:)=contourData(2,(contourBorderPosition(i)+1):(contourBorderPosition(i)+contourData_num(2,i)));
        contourBorderOneHeight=contourData_num(1,i);
        %所有点的凸包
        if (rank(contourBorderOne)==1)
            continue;
        else
            convPoints=contourBorderOne(:,convhull(contourBorderOne'));
        end
        
        %凸包图像：figure; plot(convPoints(1,:),convPoints(2,:)); hold on; plot(contourBorderOne(1,:),contourBorderOne(2,:),'.');
        %求解包含所有点的圆，也可求椭圆https://www.zhihu.com/question/268327482
        %     Point=sum(convPoints,2)./size(convPoints,2);
        %     R0=sqrt(max(sum((convPoints-Point).^2)));
        %     DrawCylinder([Point',contourBorderOneHeight-step],[Point',contourBorderOneHeight],squareR0,100,'b',1,1);
        
        %求解最小圆
        circleParameters=MinCircle(convPoints);
        parameters(i,:)=[contourBorderOneHeight,circleParameters];
        %%验证最小圆
        %     figure
        %     hold on;%保持图像zd在原图上内
        %     R = circleParameters(1);
        %     alpha=0:pi/50:2*pi;%角度容[0,2*pi]
        %     %R=2;%半径
        %     x=R*cos(alpha)+circleParameters(2);
        %     y=R*sin(alpha)+circleParameters(3);
        %     plot(x,y,'-')
        %     axis equal
        %     plot(convPoints(1,:),convPoints(2,:)); hold on; plot(contourBorderOne(1,:),contourBorderOne(2,:),'.');
    end
    %% 优化圆柱
    %阈值：判断多个圆是否为描述同一障碍物：关联圆
    parameters(all(parameters==0,2),:) = [];             %删除全0行
    radiusThreshold=2;
    %关联圆合并
    circlesNumber=size(parameters,1);
    allPairs=[];
    for i=1:circlesNumber
        %两圆心距离+较小圆半径<较大圆半径：包含关系+半径差
        partner=find(min(parameters(:,2),parameters(i,2))+sqrt(sum((parameters(:,[3,4])-parameters(i,[3,4])).^2,2))...
            <=max(parameters(:,2),parameters(i,2))...
            & abs(parameters(:,2)-parameters(i,2))<radiusThreshold);
        if isempty(partner)
            allPairs=[allPairs;[i,i]];
        else
            allPairs=[allPairs;[i*ones(length(partner),1),partner]];
        end
    end
    allPairs=sort(allPairs,2);
    allPairs=sortrows(allPairs,1);
    allPairs=unique(allPairs,'rows'); %所有的关联点对（点指的是序号，即parameters里的第n个圆柱）
    
    %寻找所有关联圆/n步关联
    records=zeros(circlesNumber,1);
    for i=1:circlesNumber
        if records(i)==0
            %第一步
            corPairsIndex=find(allPairs(:,1)==i);
            corPairs=allPairs(corPairsIndex,:);
            %寻找n步关联圆
            groupFront=[];
            groupBack=unique(corPairs);
            while length(groupFront)~=length(groupBack) || all(sort(groupFront)~=sort(groupBack))
                corPairs=[];
                for j=1:length(groupBack)
                    corPairsIndex=[find(allPairs(:,1)==groupBack(j));find(allPairs(:,2)==groupBack(j))];
                    corPairs=[corPairs;allPairs(corPairsIndex,:)];
                end
                groupFront=groupBack;
                groupBack=unique(corPairs);
            end
            %划分圆的关联组进行合并
            group=groupBack;
            records(group,:)=1;
            %取出相交圆参数
            multiCircles=parameters(group,:);
            %圆上取点
            thetaStep=100;
            theta = rand(1,thetaStep)*2*pi;
            multiCirclesPointsX=multiCircles(:,2)*cos(theta)+multiCircles(:,3);
            multiCirclesPointsY=multiCircles(:,2)*sin(theta)+multiCircles(:,4);
            multiCirclesPoints=[reshape(multiCirclesPointsX,numel(multiCirclesPointsX),1),reshape(multiCirclesPointsY,numel(multiCirclesPointsY),1)];
            convCirclesPoints=multiCirclesPoints(convhull(multiCirclesPoints),:);
            %多点画最小圆
            singleCircleParameters=MinCircle(convCirclesPoints');
            parameters(group,:)=zeros(length(group),4);
            %多个相交圆-公切圆：高度取最大
            parameters(group(1),:)=[max(multiCircles(:,1)),singleCircleParameters];
        end
    end
    parameters(all(parameters==0,2),:) = [];             %删除全0行
    parameters(:,1)=parameters(:,1)+step;                %由于等高线步长误差，导致近似圆柱的高度比真实障碍物低N米，其中N<step
    
    %% 画图优化，将 下底面 圆心二维坐标 拓展至 三维坐标
    circlesNumber=size(parameters,1);
    parameters(:,5)=zeros(circlesNumber,1);
    for i=1:circlesNumber
        %被包含
        included=find(parameters(i,2)+sqrt(sum((parameters(:,[3,4])-parameters(i,[3,4])).^2,2))<parameters(:,2));
        if isempty(included)==0
            if parameters(i,1)>max(parameters(included,1))
                parameters(i,5)=max(parameters(included,1));
            else
                parameters(i,:)=zeros(1,5);
            end
        end
    end
    parameters(all(parameters==0,2),:) = [];             %删除全0行
    csvwrite('CylinderApproximationNew.csv',parameters); %圆柱近似地形数据矩阵
end
%% 画圆柱近似后的地形图像
figure
colorData=colormap(parula);
colorScope=size(colorData,1);
heightScope=max(parameters(:,1))-min(parameters(:,1))+1;
for i=1:size(parameters,1)
    plateColorIndex=ceil((parameters(i,1)-min(parameters(:,1))+1)/heightScope*colorScope);
    plateColor=colorData(plateColorIndex,:);
    DrawCylinder(parameters(i,[3,4,5]),[parameters(i,[3,4]),parameters(i,1)],parameters(i,2),20,plateColor,1,0);
    %scatter(parameters(i,3),parameters(i,4),'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',3);
end
%% 画原地形图像
% map = map(:, 1:320);
mapColorSet=surf(map);
% set(mapColorSet,'FaceAlpha',1,'EdgeAlpha',0.6,'LineWidth',0.2)
% screenSize=get(0,'ScreenSize');
% set(gcf,'position',[screenSize(3)/5,screenSize(4)/5,1200,800]) %图像大小
axis([0 size(map,2) 0 size(map,1) 0 max(max(map))]) %坐标轴范围
view(-20,36) %给三维空间图形设置观察点的方位角az与仰角el
%[az,el]=view 返回当前的方位角az与仰角el
%xlabel('Length/5m')
%ylabel('Width/5m')
%zlabel('Height/m')
%set(gca,'xTickLabel',{'500','400','300','200','100','0'});
%axis equal      %横坐标和纵坐标轴采用等长刻度
%set(gca,'FontSize',30,'Fontname', 'Times New Roman');%设置字体
%grid on
%set(gca,'LineWidth',1.5)

hold on
%% 圆柱化之后的 地形栅格数据
gridm=size(map,1);
gridn=size(map,2);
mapCylinder=zeros(gridm,gridn);
for i=1:gridm
    for j=1:gridn
        index=find( sum( ([j,i]-parameters(:,[3,4])).^2 , 2) < parameters(:,2).^2);
        if isempty(index)==0
            mapCylinder(i,j)=max(parameters(index,1));
        end
    end
end
csvwrite('MelbourneMapCylinder.csv',mapCylinder); %圆柱近似地形数据矩阵
%% 扩大一层保护区
parametersB=parameters;
parametersB(:,2)=parametersB(:,2)+5;
parametersB(:,1)=parametersB(:,1)+5;
gridm=size(map,1);
gridn=size(map,2);
mapCylinderB=zeros(gridm,gridn);
for i=1:gridm
    for j=1:gridn
        index=find(sum( ([j,i]-parametersB(:,[3,4])).^2 , 2) < parametersB(:,2).^2);
        if isempty(index)==0
            mapCylinderB(i,j)=max(parametersB(index,1));
        end
    end
end
csvwrite('MelbourneMapCylinderB.csv',mapCylinderB); %圆柱近似地形数据矩阵
%% 画原圆柱化后的地形栅格图像
mapColorSet=surf(mapCylinder);
set(mapColorSet,'FaceAlpha',0.4,'EdgeAlpha',0,'LineWidth',0.15)

% screenSize=get(0,'ScreenSize');
% set(gcf,'position',[screenSize(3)/5,screenSize(4)/5,1200,800]) %图像大小
axis([0 size(mapCylinder,2) 0 size(mapCylinder,1) 0 max(max(mapCylinder))]) %坐标轴范围
view(-20,36) %给三维空间图形设置观察点的方位角az与仰角el
%[az,el]=view 返回当前的方位角az与仰角el
xlabel('Length/5m')
ylabel('Width/5m')
zlabel('Height/m')
%set(gca,'xTickLabel',{'500','400','300','200','100','0'});
%axis equal      %横坐标和纵坐标轴采用等长刻度
set(gca,'FontSize',30,'Fontname', 'Times New Roman');%设置字体
grid off
set(gca,'LineWidth',1.5)

ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
ax.ZAxis.Visible='off';
%% 
figure(100);
h =30;
for i=1:size(map,1)
    for j=1:size(map,2)
        E(i,j)=round(map(i,j)/h); % 10m为单位栅格的标准高度（可调）

    end
end
h = bar3(E);
% 设置颜色为单一的淡灰色
color = [0 0 1]; % RGB值表示淡灰色
axis([0 size(E,2) 0 size(E,1) 0 max(max(E))])
% 遍历每个子图并设置颜色
for k = 1:length(h)
    h(k).FaceColor = color;
end
set(gca, 'XDir', 'reverse');
set(gca, 'YDir', 'reverse');
set(h,'edgecolor','none');
grid off
ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
ax.ZAxis.Visible='off';

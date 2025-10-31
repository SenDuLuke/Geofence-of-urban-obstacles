clear all
clc
close all
addpath(genpath(pwd));
format long
if exist('MelbourneMap.csv','file')~=0
    map=xlsread('MelbourneMap.csv');
else
    %%经纬度坐标转换
    data=xlsread('澳大利亚墨尔本局部地区.csv');
    minLatitude=min(data(:,2));
    minLongitude=min(data(:,3));
    [rowm,~]=size(data);
    X=zeros(rowm,2);
    for i=1:rowm
        [X(i,1),X(i,2)]=TransLatLon(data(i,2),data(i,3),0,minLongitude);
    end
    X(:,2)=X(:,2)-min(X(:,2));
    %%转化为栅格数据，5米间隔
    gridm=ceil(max(X(:,2))/5);
    gridn=ceil(max(X(:,1))/5);
    Y=ceil(X/5);
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
    csvwrite('MelbourneMap.csv',map); %地形数据矩阵
end
%% 等高线图
step=20; %等高线步长
contourStep=step*2:step:ceil(max(max(map))/step)*step;
visible=figure;
contourData=contour(map,contourStep);  %等高线图（contourData表示等高线值及其坐标数量+若干坐标）
 % pcolor(X,Y,Z);shading interp %伪彩色图
set(visible,'visible','off');

%% 等高线处理
if exist('CylinderApproximation.csv','file')~=0
    parameters=xlsread('CylinderApproximation.csv');
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
        
        %相交圆/类同心圆合并
        threshold=20;%设置圆柱半径不超过20个格子，即100m
        [concentric,~]=find(parameters(1:i,2)<threshold & sum((parameters(1:i,[3,4])-circleParameters([2,3])).^2,2)<(max(circleParameters(1),parameters(1:i,2)).^2));
        if isempty(concentric)==0
            %取出相交圆参数
            multiCircles=parameters(concentric,:);
            %圆上取点
            thetaStep=1000;
            theta = rand(1,thetaStep)*2*pi;
            multiCirclesPointsX=multiCircles(:,2)*cos(theta)+multiCircles(:,3);
            multiCirclesPointsY=multiCircles(:,2)*sin(theta)+multiCircles(:,4);
            multiCirclesPoints=[reshape(multiCirclesPointsX,numel(multiCirclesPointsX),1),reshape(multiCirclesPointsY,numel(multiCirclesPointsY),1)];
            convCirclesPoints=multiCirclesPoints(convhull(multiCirclesPoints),:);
            %多点画最小圆
            singleCircleParameters=MinCircle(convCirclesPoints');
            parameters(concentric,:)=zeros(length(concentric),4);
            %多个相交圆-公切圆：高度取最大
            parameters(i,:)=[max(multiCircles(:,1)),singleCircleParameters];
        end

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
    parameters(all(parameters==0,2),:) = [];           %删除全0行
    parameters(:,1)=parameters(:,1)+step;              %由于等高线步长误差，导致近似圆柱的高度比真实障碍物低N米，其中N<step
    csvwrite('CylinderApproximation.csv',parameters); %圆柱近似地形数据矩阵
end
%% 画圆柱近似后的地形图像
for i=1:size(parameters,1)
    DrawCylinder([parameters(i,[3,4]),0],[parameters(i,[3,4]),parameters(i,1)],parameters(i,2),100,'b',1,1);
    alpha(0.2)
end
%% 画原地形图像
figure
mapColorSet_origin=surf(map);
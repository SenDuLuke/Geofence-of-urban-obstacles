clear all
clc
close all
addpath(genpath(pwd));
format long
if exist('MelbourneMap.csv','file')~=0
    map=xlsread('MelbourneMap.csv');
else
    %%��γ������ת��
    data=xlsread('�Ĵ�����ī�����ֲ�����.csv');
    minLatitude=min(data(:,2));
    minLongitude=min(data(:,3));
    [rowm,~]=size(data);
    X=zeros(rowm,2);
    for i=1:rowm
        [X(i,1),X(i,2)]=TransLatLon(data(i,2),data(i,3),0,minLongitude);
    end
    X(:,2)=X(:,2)-min(X(:,2));
    %%ת��Ϊդ�����ݣ�5�׼��
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
    csvwrite('MelbourneMap.csv',map); %�������ݾ���
end
%% �ȸ���ͼ
step=20; %�ȸ��߲���
contourStep=step*2:step:ceil(max(max(map))/step)*step;
visible=figure;
contourData=contour(map,contourStep);  %�ȸ���ͼ��contourData��ʾ�ȸ���ֵ������������+�������꣩
 % pcolor(X,Y,Z);shading interp %α��ɫͼ
set(visible,'visible','off');

%% �ȸ��ߴ���
if exist('CylinderApproximation.csv','file')~=0
    parameters=xlsread('CylinderApproximation.csv');
else
    contourBorder=(rem(contourData(1,:),step)==0) & (rem(contourData(2,:),1)==0);   %Ѱ�ҵȸ���ֵ���ֽ�㣩
    [~,contourBorderPosition]=find(contourBorder==1);                               %�ֽ����contourData�е�λ��
    contourData_num=contourData(:,contourBorderPosition);                           %�ֽ���Ӧ����������
    %����ǰ������Ԫ�ز�ֵ���޳���һ��û��ɸѡ���Ľ��
    for i=2:size(contourData_num,2)
        if(contourData_num(1,i)-contourData_num(1,i-1)>step)
            contourBorder(contourBorderPosition(i))=0;
        end
    end
    %���»�ȡһ��
    [~,contourBorderPosition]=find(contourBorder==1);
    contourData_num=contourData(:,contourBorderPosition);
    
    %%��ȡ�����ȸ�������
    parameters=zeros(size(contourData_num,2),4); %�������ܣ��߶ȡ��뾶���������꣩
    for i=1:size(contourData_num,2)
        contourBorderOne=zeros(2,contourData_num(2,i));
        contourBorderOne(1,:)=contourData(1,(contourBorderPosition(i)+1):(contourBorderPosition(i)+contourData_num(2,i)));
        contourBorderOne(2,:)=contourData(2,(contourBorderPosition(i)+1):(contourBorderPosition(i)+contourData_num(2,i)));
        contourBorderOneHeight=contourData_num(1,i);
        %���е��͹��
        if (rank(contourBorderOne)==1)
            continue;
        else
            convPoints=contourBorderOne(:,convhull(contourBorderOne'));
        end
        
        %͹��ͼ��figure; plot(convPoints(1,:),convPoints(2,:)); hold on; plot(contourBorderOne(1,:),contourBorderOne(2,:),'.');
        %���������е��Բ��Ҳ������Բhttps://www.zhihu.com/question/268327482
        %     Point=sum(convPoints,2)./size(convPoints,2);
        %     R0=sqrt(max(sum((convPoints-Point).^2)));
        %     DrawCylinder([Point',contourBorderOneHeight-step],[Point',contourBorderOneHeight],squareR0,100,'b',1,1);
        
        %�����СԲ
        circleParameters=MinCircle(convPoints);
        parameters(i,:)=[contourBorderOneHeight,circleParameters];
        
        %�ཻԲ/��ͬ��Բ�ϲ�
        threshold=20;%����Բ���뾶������20�����ӣ���100m
        [concentric,~]=find(parameters(1:i,2)<threshold & sum((parameters(1:i,[3,4])-circleParameters([2,3])).^2,2)<(max(circleParameters(1),parameters(1:i,2)).^2));
        if isempty(concentric)==0
            %ȡ���ཻԲ����
            multiCircles=parameters(concentric,:);
            %Բ��ȡ��
            thetaStep=1000;
            theta = rand(1,thetaStep)*2*pi;
            multiCirclesPointsX=multiCircles(:,2)*cos(theta)+multiCircles(:,3);
            multiCirclesPointsY=multiCircles(:,2)*sin(theta)+multiCircles(:,4);
            multiCirclesPoints=[reshape(multiCirclesPointsX,numel(multiCirclesPointsX),1),reshape(multiCirclesPointsY,numel(multiCirclesPointsY),1)];
            convCirclesPoints=multiCirclesPoints(convhull(multiCirclesPoints),:);
            %��㻭��СԲ
            singleCircleParameters=MinCircle(convCirclesPoints');
            parameters(concentric,:)=zeros(length(concentric),4);
            %����ཻԲ-����Բ���߶�ȡ���
            parameters(i,:)=[max(multiCircles(:,1)),singleCircleParameters];
        end

        %%��֤��СԲ
        %     figure
        %     hold on;%����ͼ��zd��ԭͼ����
        %     R = circleParameters(1);
        %     alpha=0:pi/50:2*pi;%�Ƕ���[0,2*pi]
        %     %R=2;%�뾶
        %     x=R*cos(alpha)+circleParameters(2);
        %     y=R*sin(alpha)+circleParameters(3);
        %     plot(x,y,'-')
        %     axis equal
        %     plot(convPoints(1,:),convPoints(2,:)); hold on; plot(contourBorderOne(1,:),contourBorderOne(2,:),'.');
    end    
    parameters(all(parameters==0,2),:) = [];           %ɾ��ȫ0��
    parameters(:,1)=parameters(:,1)+step;              %���ڵȸ��߲��������½���Բ���ĸ߶ȱ���ʵ�ϰ����N�ף�����N<step
    csvwrite('CylinderApproximation.csv',parameters); %Բ�����Ƶ������ݾ���
end
%% ��Բ�����ƺ�ĵ���ͼ��
for i=1:size(parameters,1)
    DrawCylinder([parameters(i,[3,4]),0],[parameters(i,[3,4]),parameters(i,1)],parameters(i,2),100,'b',1,1);
    alpha(0.2)
end
%% ��ԭ����ͼ��
figure
mapColorSet_origin=surf(map);
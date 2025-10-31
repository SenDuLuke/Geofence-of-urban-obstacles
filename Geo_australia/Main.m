clear
clc
close all
addpath(genpath(pwd));
format long

%% �����ͼ
if exist('MelbourneMap.csv','file')~=0
    map=xlsread('MelbourneMap.csv');
else
% ��γ������ת��
    data=xlsread('�Ĵ�����ī�����ֲ�����.csv');
    minLatitude=min(data(:,2));
    minLongitude=min(data(:,3));
    [rowm,~]=size(data);
    X=zeros(rowm,2);
    for i=1:rowm
        [X(i,1),X(i,2)]=TransLatLon(data(i,2),data(i,3),0,minLongitude);
    end
%%     
% ת��Ϊդ�����ݣ�5�׼��
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
    % csvwrite('MelbourneMap.csv',map); %�������ݾ���
end

%% �ȸ���ͼ
step=10; %�ȸ��߲���
contourStep=step:step:ceil(max(max(map))/step)*step;
visible=figure;
contourData=contour(map,contourStep);  %�ȸ���ͼ��contourData��ʾ�ȸ���ֵ������������+�������꣩
%pcolor(X,Y,Z);shading interp %α��ɫͼ
%set(visible,'visible','off');

%% �ȸ��ߴ���
if exist('CylinderApproximationNew.csv','file')~=0
    parameters=xlsread('CylinderApproximationNew.csv');
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
    %% �Ż�Բ��
    %��ֵ���ж϶��Բ�Ƿ�Ϊ����ͬһ�ϰ������Բ
    parameters(all(parameters==0,2),:) = [];             %ɾ��ȫ0��
    radiusThreshold=2;
    %����Բ�ϲ�
    circlesNumber=size(parameters,1);
    allPairs=[];
    for i=1:circlesNumber
        %��Բ�ľ���+��СԲ�뾶<�ϴ�Բ�뾶��������ϵ+�뾶��
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
    allPairs=unique(allPairs,'rows'); %���еĹ�����ԣ���ָ������ţ���parameters��ĵ�n��Բ����
    
    %Ѱ�����й���Բ/n������
    records=zeros(circlesNumber,1);
    for i=1:circlesNumber
        if records(i)==0
            %��һ��
            corPairsIndex=find(allPairs(:,1)==i);
            corPairs=allPairs(corPairsIndex,:);
            %Ѱ��n������Բ
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
            %����Բ�Ĺ�������кϲ�
            group=groupBack;
            records(group,:)=1;
            %ȡ���ཻԲ����
            multiCircles=parameters(group,:);
            %Բ��ȡ��
            thetaStep=100;
            theta = rand(1,thetaStep)*2*pi;
            multiCirclesPointsX=multiCircles(:,2)*cos(theta)+multiCircles(:,3);
            multiCirclesPointsY=multiCircles(:,2)*sin(theta)+multiCircles(:,4);
            multiCirclesPoints=[reshape(multiCirclesPointsX,numel(multiCirclesPointsX),1),reshape(multiCirclesPointsY,numel(multiCirclesPointsY),1)];
            convCirclesPoints=multiCirclesPoints(convhull(multiCirclesPoints),:);
            %��㻭��СԲ
            singleCircleParameters=MinCircle(convCirclesPoints');
            parameters(group,:)=zeros(length(group),4);
            %����ཻԲ-����Բ���߶�ȡ���
            parameters(group(1),:)=[max(multiCircles(:,1)),singleCircleParameters];
        end
    end
    parameters(all(parameters==0,2),:) = [];             %ɾ��ȫ0��
    parameters(:,1)=parameters(:,1)+step;                %���ڵȸ��߲��������½���Բ���ĸ߶ȱ���ʵ�ϰ����N�ף�����N<step
    
    %% ��ͼ�Ż����� �µ��� Բ�Ķ�ά���� ��չ�� ��ά����
    circlesNumber=size(parameters,1);
    parameters(:,5)=zeros(circlesNumber,1);
    for i=1:circlesNumber
        %������
        included=find(parameters(i,2)+sqrt(sum((parameters(:,[3,4])-parameters(i,[3,4])).^2,2))<parameters(:,2));
        if isempty(included)==0
            if parameters(i,1)>max(parameters(included,1))
                parameters(i,5)=max(parameters(included,1));
            else
                parameters(i,:)=zeros(1,5);
            end
        end
    end
    parameters(all(parameters==0,2),:) = [];             %ɾ��ȫ0��
    csvwrite('CylinderApproximationNew.csv',parameters); %Բ�����Ƶ������ݾ���
end
%% ��Բ�����ƺ�ĵ���ͼ��
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
%% ��ԭ����ͼ��
% map = map(:, 1:320);
mapColorSet=surf(map);
% set(mapColorSet,'FaceAlpha',1,'EdgeAlpha',0.6,'LineWidth',0.2)
% screenSize=get(0,'ScreenSize');
% set(gcf,'position',[screenSize(3)/5,screenSize(4)/5,1200,800]) %ͼ���С
axis([0 size(map,2) 0 size(map,1) 0 max(max(map))]) %�����᷶Χ
view(-20,36) %����ά�ռ�ͼ�����ù۲��ķ�λ��az������el
%[az,el]=view ���ص�ǰ�ķ�λ��az������el
%xlabel('Length/5m')
%ylabel('Width/5m')
%zlabel('Height/m')
%set(gca,'xTickLabel',{'500','400','300','200','100','0'});
%axis equal      %�����������������õȳ��̶�
%set(gca,'FontSize',30,'Fontname', 'Times New Roman');%��������
%grid on
%set(gca,'LineWidth',1.5)

hold on
%% Բ����֮��� ����դ������
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
csvwrite('MelbourneMapCylinder.csv',mapCylinder); %Բ�����Ƶ������ݾ���
%% ����һ�㱣����
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
csvwrite('MelbourneMapCylinderB.csv',mapCylinderB); %Բ�����Ƶ������ݾ���
%% ��ԭԲ������ĵ���դ��ͼ��
mapColorSet=surf(mapCylinder);
set(mapColorSet,'FaceAlpha',0.4,'EdgeAlpha',0,'LineWidth',0.15)

% screenSize=get(0,'ScreenSize');
% set(gcf,'position',[screenSize(3)/5,screenSize(4)/5,1200,800]) %ͼ���С
axis([0 size(mapCylinder,2) 0 size(mapCylinder,1) 0 max(max(mapCylinder))]) %�����᷶Χ
view(-20,36) %����ά�ռ�ͼ�����ù۲��ķ�λ��az������el
%[az,el]=view ���ص�ǰ�ķ�λ��az������el
xlabel('Length/5m')
ylabel('Width/5m')
zlabel('Height/m')
%set(gca,'xTickLabel',{'500','400','300','200','100','0'});
%axis equal      %�����������������õȳ��̶�
set(gca,'FontSize',30,'Fontname', 'Times New Roman');%��������
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
        E(i,j)=round(map(i,j)/h); % 10mΪ��λդ��ı�׼�߶ȣ��ɵ���

    end
end
h = bar3(E);
% ������ɫΪ��һ�ĵ���ɫ
color = [0 0 1]; % RGBֵ��ʾ����ɫ
axis([0 size(E,2) 0 size(E,1) 0 max(max(E))])
% ����ÿ����ͼ��������ɫ
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

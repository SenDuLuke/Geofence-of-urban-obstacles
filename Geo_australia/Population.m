clear all
clc
close all
addpath(genpath(pwd));
format long

%% ���� �߶�դ�� �� ¥���ĵ� ����
map=xlsread('MelbourneMap.csv');
cy_map=xlsread('MelbourneMapCylinder.csv');
cylinder=xlsread('CylinderApproximationNew.csv');

%  %����¥���ĵ�λ��
%   mesh(cy_map);
%    hold on 
%   for i=1:length(cylinder(:,1))
%    scatter(cylinder(i,3),cylinder(i,4));
%     hold on
%   end

%% �˿��ܶȲ�������
p_ave=0.8;

% ���㽨�����Բ����� Ϊ�˿ڲ�����׼
for i=1:length(cylinder(:,1))
    if  cylinder(i,2)<5
        cylinder(i,6)=0;
        
    else
        cylinder(i,6)=(cylinder(i,1)-cylinder(i,5));
    end
end

for i=1:length(cylinder(:,1))
    
    cylinder(i,7)= p_ave*cylinder(i,6);
    
    
end

%% �˿��ܶ�դ�����
 gridm=size(map,1);
 gridn=size(map,2);
 num_build=length(cylinder(:,1));
 population_Mel=zeros(gridm,gridn);
 
%  p_index=zeros(length(cylinder(:,1)),1);
 
 for k=1:num_build
     
     
   for i=1:gridm
       for j=1:gridn
           
              mxtr=5 * ([j,i]-[cylinder(k,3),cylinder(k,4)]);
              dis= sqrt(mxtr*mxtr')/1000;
              
          if dis<0.2
              
                 population_Mel(i,j)= cylinder(k,7) * exp(1-dis^2);
                 
          end
          
       end
          
          
              
              
   end
   
  end           
          
%% 
  mesh (population_Mel)
% %% 
% y=zeros(50,1);
% x=zeros(50,1);
% 
% for i=1:50
%     y(i)= 100*exp(1-i^2);
%     x(i)=i;
% end
% 
% plot(x,y);  
          
          
          
          
          
          
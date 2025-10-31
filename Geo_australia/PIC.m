clear all
clc
close all
addpath(genpath(pwd));
format long
X=xlsread('Result_Data_2.xls','1000');
[m,n]=size(X);
t0=1/20;
t=t0/2:t0:m*t0;
figure
plot(t,X(:,5),'b')   %低空空域安全态势时间曲线图
figure
x=t;
y=X(:,5).';
sizeMarker = linspace(1, 100, length(x));    % 比0大，值越大标记越大
colorMarker = y;   % 颜色渐变
% subplot(1,2,1)
% scatter(x, y, sizeMarker, colorMarker, 'o', 'filled')
% subplot(1,2,2)
patch([x NaN],[y NaN],[colorMarker NaN],'Marker','o','EdgeColor','interp','MarkerFaceColor','flat')
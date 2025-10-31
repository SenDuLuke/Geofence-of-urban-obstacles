function [X,Y]=TransLatLon(lat,lon,B0,L0)
% 调用格式：[X,Y]=TransLatLon(lat,lon,B0,L0)
% 输入： lat为纬度，lon为经度（单位为弧度），B0：坐标原点的纬度，对于地球必须设置为0，L0：坐标原点的经度，单位为弧度((rad);
% 输出：X为水平直角坐标， Y为纵向直角坐标，单位为米(m);
B=lat*pi/180;L=lon*pi/180;
B0=B0*pi/180;L0=L0*pi/180;%弧度转换为度数
a=6378137; b=6356752.3142;% a:地球椭球体长半轴，b:地球椭球体短半轴，单位为米(m);
e=sqrt(1-(b/a)^2);ee=sqrt((a/b)^2-1);%偏心率
K=(((a^2*cos(B0))/b)/sqrt(1+(ee)^2*(cos(B0))^2));
X=K*(L-L0); 
Y=K*log(tan(pi/4+B/2)*((1-e*sin(B))/(1+e*sin(B)))^(e/2));
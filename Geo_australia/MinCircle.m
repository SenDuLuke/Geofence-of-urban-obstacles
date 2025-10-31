function circleParameters = MinCircle( convPoints )
detSigma=@(x)x(1)^2; %三个参数分别为半径、中心点坐标
c = @(x)sum((convPoints'-[x(2),x(3)]).^2,2)-x(1)^2;
ceq = @(x)[];
nonlinfcn = @(x)deal(c(x),ceq(x));
%初值半径
initialPoint=sum(convPoints,2)./size(convPoints,2);
initialR0=sqrt(max(sum((convPoints-initialPoint).^2)));
%求解
circleParameters=fmincon(detSigma,[initialR0,initialPoint'],[],[],[],[],[0,0,0],[],nonlinfcn);
end


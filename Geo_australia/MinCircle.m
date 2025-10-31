function circleParameters = MinCircle( convPoints )
detSigma=@(x)x(1)^2; %���������ֱ�Ϊ�뾶�����ĵ�����
c = @(x)sum((convPoints'-[x(2),x(3)]).^2,2)-x(1)^2;
ceq = @(x)[];
nonlinfcn = @(x)deal(c(x),ceq(x));
%��ֵ�뾶
initialPoint=sum(convPoints,2)./size(convPoints,2);
initialR0=sqrt(max(sum((convPoints-initialPoint).^2)));
%���
circleParameters=fmincon(detSigma,[initialR0,initialPoint'],[],[],[],[],[0,0,0],[],nonlinfcn);
end


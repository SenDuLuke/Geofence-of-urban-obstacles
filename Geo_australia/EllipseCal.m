    %���������е��Բ��Ҳ������Բhttps://www.zhihu.com/question/268327482
    detSigma=@(x)x(1)^2;
    c = @(x)(contourBorderOne'-[x(2),x(3)])*[x(1),0;0,x(1)]^(-1)*(contourBorderOne'-[x(2),x(3)])'-1;
    ceq = @(x)[];
    nonlinfcn = @(x)deal(c(x),ceq(x));
    %��ֵ�뾶
    initialPoint=sum(contourBorderOne,2)./size(contourBorderOne,2);
    initialR0=1/max(sum((contourBorderOne-initialPoint).^2));
    %���
    CircleParameters=fmincon(detSigma,[initialR0,initialPoint'],[],[],[],[],[],[],nonlinfcn);
function [x,minf]=minMixFun(f,g,h,x0,r0,c,var,eps)
while 1
    FF=r0*FE+FH/sqrt(r0);
    SumF=f+FF;
    [x2,minf]=minNT(SumF,transpose(x1),var);
    if norm(x2-x1)<=eps
        x=x2;
        break;
    else
        r0=c*r0;
        x1=x2;
    end
end
minf=Funval(f,var,x);
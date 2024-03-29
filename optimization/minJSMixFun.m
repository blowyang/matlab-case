function [x,minf]=minJSMixFun(f,g,h,x0,r0,c,var,eps)
while 1
    FF=r0*FE+FH/sqrt(r0);
    SumF=f+FF;
    a0=(c*x1-x2)/(c-1);
    x2=a0+(x1-a0)*c^2;
    [x3,minf]=minNT(SumF,transpose(x2),var);
    if norm(x3-x2)<=eps
        x=x3;
        break;
    else
        r0=c*r0;
        x1=x2;
        x2=x3;
    end
end
minf=Funval(f,var,x);
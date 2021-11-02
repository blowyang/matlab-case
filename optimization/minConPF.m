function[x,minf]=minConPF(f,x0,g,h,c1,p,var,eps)
format long;
if nargin ==7
    eps=1.0e-5;
end
k=0;
FE=0;
for i=1:length(h)
	FE=FE+(h(i))^2;
end
x1=transpose(x0);
x2=inf;
while 1
    M=c1*p;
	FF=M*FE;
	gx= Funval(g,var,x1);
	gF=0;
	for i=1:length(h)
        if gx(i)<0
        gF=gF+M*(g(i)^2);
        end
    end
	SumF=f+FF+gF;
	[x2,minf]=minNT(SumF,transpose(x1),var);
	if norm(x2-x1)<=eps
        x=x2;
        break;
    else
        c1=M;
        x1=x2;
	end
end
minf= Funval(f,var,x);
format short
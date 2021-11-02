function [x,minf]=minNT(f,x0,var,eps)
format long;
if nargin==3
	eps=1.0e-5;
end
tol=1;
x0=transpose(x0);
while tol>eps
	gradf=jacobian(f,var);
	jacf=jacobian(gradf,var);
	v =Funval(gradf,var,x0);
	tol= norm(v);
	pv= Funval(jacf,var,x0);
	p=-inv(pv)*transpose(v);
	x1=x0+p;
	x0=x1;
end
x=x1;
minf=Funval(f,var,x);
format short
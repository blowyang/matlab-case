 clear
 clc
syms s t;
f= s^2-s*t+t-s+1;
g=[s^2+t^2-6;s;t];h=[2*s+3*t-9];
[x,minf]=minConPF(f,[2 2],g,h,0.05,2,[s t])
syms x1 x2;
tic
f=-x1+x2;
g=[log(x2);x1;x2];
h=[x1+x2-1];
[x,minf]=minMixFun(f,g,h,[2 2],10,0.5,[x1,x2],1.0e-5);
[x,minf]=minJSMixFun(f,g,h,[2 2],10,0.5,[x1,x2],1.0e-5);
toc
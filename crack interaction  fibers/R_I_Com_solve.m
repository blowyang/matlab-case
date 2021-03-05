%形成复系数的方程后，将实部和虚部分开处理,该程序仅能处理，Z1、Z2、Zd位于坐标轴x上的情况
clear;
clc;

S11=0;
S22=0;
S12=0;
E3=60;
%E3=0;
B1=(S11+S22)/4;
C1=0;
B2=(S22-S11)/2;
C2=S12;

Z1=-50;
% Zd=2.2;
Z2=50;
R1=1;
R2=1;

bx=0;
by=1;

e31=-6.5;

vm=0.3;
vf=0.3;
km=3-4*vm;
kf=3-4*vf;

Ef=12.6e6;
Em=2*Ef;

um=Em/2/(1+vm);
uf=Ef/2/(1+vf);

CE=um/uf;

gama=um*(by-1i*bx)/pi/(1+km);
da=0.2;
st=-4;
J=39;%计算迭代次数
N=50;%迭代次数


for kd=1:1:J
    Zd=st+kd*da;
    [A,B]=xishu(N,Z1,Z2,Zd,R1,R2,S11,S12,S22,E3,bx,by,um,CE);
    [CA, CB]=Cxishu(A ,B,N);
    
    CB=CB.';
    CX=CA\CB;
    for k=1:1:8*N
        tem_1=CX(2*k-1);
        X(k)=tem_1;
    end
    
    RX=real(X);
    IX=imag(X);
   
    Force_1(kd)=FX(N,R1,R2,Z1,Z2,RX,IX,Zd,bx,by)*pi*R2*(1+km)/(um*bx*bx);
    Force_2(kd)=FY(N,R1,R2,Z1,Z2,RX,IX,Zd,bx,by)*pi*R2*(1+km)/(um*bx*bx);
end
% clear A;    
% Clear B;
% clear CA;
% Clear CB;    
% clear CX;
% clear X;
% clear RX;
% clear IX;
Force_1=Force_1';
Force_2=Force_2';


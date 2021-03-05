%形成复系数的方程后，将实部和虚部分开处理,该程序仅能处理，Z1、Z2、Zd位于坐标轴x上的情况
clear;
clc;

S11=0;
e31=-6.5;
S12=0;

%E3=0;
Z1=-2;
Z2=2;
R1=1;
R2=1;

bx=0;

vm=0.3;
vf=0.3;
km=3-4*vm;
kf=3-4*vf;

Ef=12.6e10;
Em=0.1*Ef;

um=Em/2/(1+vm);
uf=Ef/2/(1+vf);
% 
% CE=um/uf;
%
N1=50; 
N2=50;
N=50;
d=0;
c=0.25;

S22=0;
E3=0;
by=1;

gama=um*(by-1i*bx)/pi/(1+km);


for kd=1:1:N    
    Zd=c*cos((2*kd-1)*pi/2/N)+d; 
    [A,B]=xishu(N2,Z1,Z2,Zd,R1,R2,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef);
    i=1;
    for n=1:1:8*N2
        j=1;
        m=1;
        for k=1:1:8*N2
            A1(i,m)=A(i,j);
            j=j+1;
            A2(i,m)=A(i,j);
            j=j+1;
            m=m+1;
        end
        i=i+1;
    end
    RA=A1+A2;
    IA=A1-A2;
    B=B.';
    RB=real(B);
    IB=imag(B);   
    RX=RA\RB;
    IX=IA\IB;
    for r=1:1:N-1
        Z=c*cos(r*pi/N)+d;         
        AA(r,kd)=c*pi/N*Str22(S11,S22,S12,N2,R1,R2,Z1,Z2,RX,IX,Z,Zd,gama);
    end    
end
for kd=1:1:N
    AA(N,kd)=1;
end

by=0;
S22=1e6;
CS22=S22;
E3=-0.1*S22/e31;
gama=um*(by-1i*bx)/pi/(1+km);

for kd=-10:1:10
    E3=kd*(-0.5*S22/e31);
    [A,B]=xishu(N1,Z1,Z2,Zd,R1,R2,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef);
    i=1;
    for n=1:1:8*N1
        j=1;
        m=1;
        for k=1:1:8*N1
            A1(i,m)=A(i,j);
            j=j+1;
            A2(i,m)=A(i,j);
            j=j+1;
            m=m+1;
        end
        i=i+1;
    end
    RA=A1+A2;
    IA=A1-A2;
    B=B.';
    RB=real(B);
    IB=imag(B);
    RX=RA\RB;
    IX=IA\IB;
    %%%%%%%%%%%
    
    for r=1:1:N-1
        Z=c*cos(r*pi/N)+d;
        P(r)=Z;
        BB(r)=-Str22(S11,S22,S12,N1,R1,R2,Z1,Z2,RX,IX,Z,Zd,gama);
    end
    BB(N)=0;
    BB=BB.';
    F=AA\BB;
    BB=BB';
    KI(kd+11)=-2*um*F(N)/(1+km)/CS22;
end


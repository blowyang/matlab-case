%形成复系数的方程后，将实部和虚部分开处理,该程序仅能处理，Z1、Z2、Zd位于坐标轴x上的情况
clear;
clc;

S11=0;
e31=-6.5;
S12=0;
R1=1;
R2=1;
R3=1;

d=0;%裂纹中心所在坐标
c=0.25;%半裂纹的长度

vm=0.3;
vf=0.3;
km=3-4*vm;
kf=3-4*vf;

Ef=12.6e10;
Em=10*Ef;

um=Em/2/(1+vm);
uf=Ef/2/(1+vf);
% 
% CE=um/uf;
%

N1=20; 
N2=20;
N=50;
da=0.2;
st=-0.25;
for count=1:1:1
    Zd=20;%位所在位置坐标错
    by=0;
    bx=0;
    S22=1e6;
    CS22=S22;
    E3=2*S22/e31;
    gama=um*(by-1i*bx)/pi/(1+km);
%     Z1=-st-da*count;
%     Z2=st+da*count;
%     Z3=100;
    
    Z1=st-1.5+da*count-1.5*1i;
    Z2=st+1.5+da*count;
    Z3=st-1.5+da*count+1.5*1i;
    [A,B]=xishu(N1,Z1,Z2,Z3,Zd,R1,R2,R3,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef);
    [CA, CB]=Cxishu(A ,B,N1);
    clear A;
    clear B;
    CB=CB.';
    CX=CA\CB;
    clear CA;
    clear CB;
    X=get_xishu(CX,N1);
    clear CX;
    
    %%%%%%%%%%%
    
    for r=1:1:N-1       
        Z=c*cos(r*pi/N)+d;
        BB(r)=-Str22(S11,S22,S12,N1,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama);
%         BB(N+r)=-Str12(S11,S22,S12,N1,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama);
    end
    BB(N)=0;
%     BB(2*N)=0;
    clear RX;
    clear IX;
    
    S22=0;
    E3=0;
    for kd=1:1:N  
        disp('count');
        disp(count);
        disp('kd');
        disp(kd);
        bx=0;
        by=1;
        gama1=um*(by-1i*bx)/pi/(1+km);   
        Zd=c*cos((2*kd-1)*pi/2/N)+d; 
        [A,B]=xishu(N2,Z1,Z2,Z3,Zd,R1,R2,R3,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef);
        [CA, CB]=Cxishu(A ,B,N2);
        clear A;
        clear B;
        
        CB=CB.';
        CX=CA\CB;
        clear CA;
        clear CB;
        X1=get_xishu(CX,N2);
        clear CX;
        
%         by=0;
%         bx=1;
%         gama2=um*(by-1i*bx)/pi/(1+km);   
%         [A,B]=xishu(N2,Z1,Z2,Z3,Zd,R1,R2,R3,S11,S12,S22,E3,bx,by,e31,vm,vf,Em,Ef);
%         [CA, CB]=Cxishu(A ,B,N2);
%         clear A;
%         clear B;
%         
%         CB=CB.';
%         CX=CA\CB;
%         clear CA;
%         clear CB;
%         X2=get_xishu(CX,N2);
%         clear CX;
% 
        for r=1:1:N-1
            Z=c*cos(r*pi/N)+d; 
            AA(r,kd)=c*pi/N*Str22(S11,S22,S12,N1,R1,R2,R3,Z1,Z2,Z3,X1,Z,Zd,gama1);
%             AA(r,N+kd)=c*pi/N*Str22(S11,S22,S12,N2,R1,R2,R3,Z1,Z2,Z3,X2,Z,Zd,gama2);
%             AA(N+r,kd)=c*pi/N*Str12(S11,S22,S12,N1,R1,R2,R3,Z1,Z2,Z3,X1,Z,Zd,gama1);
%             AA(N+r,N+kd)=c*pi/N*Str12(S11,S22,S12,N2,R1,R2,R3,Z1,Z2,Z3,X2,Z,Zd,gama2);
        end
    end
    for kd=1:1:N
        AA(N,kd)=1;
%         AA(N,N+kd)=0;
%         AA(2*N,kd)=0;
%         AA(2*N,N+kd)=1;
    end
    BB=BB.';
    F=AA\BB;
    BB=BB.';
    KI_left(count)=-2*um*(F(N))/(1+km)/CS22;
    KI_right(count)=2*um*(F(1))/(1+km)/CS22;
%     KII_left(count)=-2*um*(F(2*N))/(1+km)/CS22;
%     KII_right(count)=2*um*(F(N+1))/(1+km)/CS22;
  
end


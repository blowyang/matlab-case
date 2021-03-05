clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
N=20;
M=15;
K=80;

C11E=12.6;
C12E=5.5;
Ef=C11E-C12E/C11E;
muf=C12E/C11E;

Gf=Ef/2/(1+muf);
kf=3-4*muf;
Gm=5*Gf/10;
km=1.8;
a=4;

e31=-6.5;
%E3infinity=0.2;
E3infinity=1;
S022=0;
S011=0;
S012=0;
% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
Em=2.6*Gm;


% vf=0.5;
% R=a*(vf/pi)^0.5;
% [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km);
% B_Matrix=Get_B_Matrix(N,M,km,R,a);
% ELoad=Get_ELoad(N,e31,E3infinity,R);
% [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm);
% H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
% [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,C_Matrix,F_Matrix,ELoad);
% 
% [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
%     Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);

NP=40;
del=0.02;
Stif=zeros(NP,2);
for np=1:NP
    vf=del*np;
    R=a*(vf/pi)^0.5;
    
    [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km);
    B_Matrix=Get_B_Matrix(N,M,km,R,a);
    ELoad=Get_ELoad(N,e31,E3infinity,R);
    [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm);
    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
    [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,C_Matrix,F_Matrix,ELoad);
    
    [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
        Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
    Stif(np,1)=vf;
    
    Stif(np,2)=delta3/a;
   
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
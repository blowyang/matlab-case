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
K=60;

C11E=12.6;
C12E=7.95;
Ef=C11E-C12E/C11E;
muf=C12E/C11E;

Gf=Ef/2/(1+muf);
kf=3-4*muf;
Gm=10*Gf;
km=1.8;
a=4;

e31=-6.5;
%E3infinity=0.2;
E3infinity=0;
S022=0;
S011=1;
S012=0;
% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
Em=2.6*Gm;

NP=20;
del=0.02;
Stif=zeros(NP,8);
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
    
    %S11,S22
    tem=delta_cof(1,2);
    Stif(np,2)=tem*a/Em;
    %S12
    tem=delta_cof(1,3);
    Stif(np,3)=tem*a/Em;
    %S33
    tem=delta_cof(2,1);
    Stif(np,4)=2*tem*a/Gm;
    %S12,S21,S13,S23
    tem=delta_cof(2,2);
    Stif(np,5)=tem*a/Gm;
    tem=delta_cof(2,3);
    Stif(np,6)=tem*a/Gm;
    tem=delta_cof(4,2);
    Stif(np,7)=tem*a/Gm;
    tem=delta_cof(4,3);
    Stif(np,8)=tem*a/Gm;
end
% for np=1:NP
%     z=R+(np-1)*del;
%     [Sm22,~,~]=Get_Matrix_Stress(a,R,Xce,Xdf,z);
%     S22(np,1)=z;
%     S22(np,2)=Sm22;
% end
% 
% S12=zeros(NP,3);
% for np=1:NP
%     zk_AB=np*a/(NP+1)-a*(1-1i)/2;
%     zk_DB=np*a*1i/(NP+1)+a*(1-1i)/2;
%     [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce,Xdf,zk_AB);
%     S12(np,1)=-a/2+np*a/(NP+1);
%     S12(np,2)=Sm12;
%     [~,~,Sm12]=Get_Matrix_Stress(a,R,Xce,Xdf,zk_DB);
%     S12(np,3)=Sm12;
% end
%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
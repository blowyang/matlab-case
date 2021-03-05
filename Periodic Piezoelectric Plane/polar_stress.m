clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
N=25;
M=20;
K=100;
Gm=10;
Gf=1;
kf=1.8;
km=1.8;
R=1;

e31=-6.5;
%E3infinity=0.2;
E3infinity=0.03;
S022=1;
S011=0;
S012=0;
vf=0.02;
    
a=R*(pi/vf)^0.5;
[E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km);
B_Matrix=Get_B_Matrix(N,M,km,R,a);
ELoad=Get_ELoad(N,e31,E3infinity,R);
[C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
[Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,C_Matrix,F_Matrix,ELoad);

[delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
    Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
Xce=Xce0+delta1*Xce1+delta2*Xce2+delta3*Xce3+delta4*Xce4;
Xdf=Xdf0+delta1*Xdf1+delta2*Xdf2+delta3*Xdf3+delta4*Xdf4;


NP=100;
del=2*pi/(NP+1);
stress=zeros(NP,3);
for np=1:NP
    z=R*exp(1i*np*del);
    
    [Smss,Smrr,~]=Get_Matrix_Polar_Stress(a,R,Xce,Xdf,z);
    stress(np,1)=np*del;
    stress(np,2)=Smss;
    stress(np,3)=Smrr;
%  
%     S13_DB_Sub_CA(np,1)=S13_DB(np,1)-S13_CA(np,1);
%     S23_AB_Sub_CD(np,1)=S23_AB(np,1)-S23_CD(np,1);
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
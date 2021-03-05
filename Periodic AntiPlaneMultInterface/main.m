clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
R2=1;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=1250;
gama1=0.1;
G1=G0/gama1;
gama2=0.1;
G2=G1/gama2;
R1=R2+0.001;
lambda=R2/R1;
lambda1=R1/a ;

N=30;
M=25;
K=20;
[E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
B_Matrix=Get_B_Matrix(N,M,gama1,lambda1);

C_Matrix=Get_C_Matrix(K,N,a,R1);
D_Matrix=Get_D_Matrix(K,M,a);

H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

[Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
S013_Divided_G0=0.1;
S023_Divided_G0=0.1;
[delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
delta1=delta1_Divided_a*a;
delta2=delta2_Divided_a*a;
Xa=delta1*Xa2+delta2*Xa1;
Xd=delta1*Xd2+delta2*Xd1;
Xe=delta1*Xe2+delta2*Xe1;

NP=100;
del=2*pi/(NP+1);
SigmaMatrix_zs=zeros(NP,1);
SigmaMatrix_zr=zeros(NP,1);
SigmaFiber_zs=zeros(NP,1);
SigmaFiber_zr=zeros(NP,1);
SigmaFiber_23=zeros(NP,1);
SigmaFiber_13=zeros(NP,1);
for np=1:NP
    zM_np=R1*cos(del*np)+R1*sin(del*np)*1i;
    zF_np=R2*cos(del*np)+R2*sin(del*np)*1i;
    [Sigma_zs_Divided_G0,Sigma_zr_Divided_G0]=Get_SigmaMatrix_Divided_G0(Xd,Xe,a,R1,zM_np);
    SigmaMatrix_zs(np,1)=Sigma_zs_Divided_G0;
    SigmaMatrix_zr(np,1)=Sigma_zr_Divided_G0;
    %返回的应力都是除以G2后的值
    [Sigma_zs_Divided_G2,Sigma_zr_Divided_G2,Sigma_23_Divided_G2,Sigma_13_Divided_G2]=Get_SigmaFiber_Divided_G2(Xa,R2,zF_np);
    %返回的应力都是除以G0后的值
    SigmaFiber_zs(np,1)=Sigma_zs_Divided_G2/gama1/gama2;
    SigmaFiber_zr(np,1)=Sigma_zr_Divided_G2/gama1/gama2; 
    SigmaFiber_23(np,1)=Sigma_23_Divided_G2/gama1/gama2;
    SigmaFiber_13(np,1)=Sigma_13_Divided_G2/gama1/gama2;
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

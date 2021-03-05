clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
%在固定G0和G2的比例情况下,研究G1变化对复合材料综合性能的影响

R2=1;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=100;
G2=1;
R1=R2+0.2;
lambda=R2/R1;
lambda1=R1/a;

N=30;
M=25;
K=20;


NP=20;
del=6/(NP+1);
log_G1DG0=zeros(NP,1);
GeDG0=zeros(NP,1);
for np=1:NP
    G1=G0*power(10,-3+del*np);
    gama1=G0/G1;
    gama2=G1/G2;
    [E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
    B_Matrix=Get_B_Matrix(N,M,gama1,lambda1);

    C_Matrix=Get_C_Matrix(K,N,a,R1);
    D_Matrix=Get_D_Matrix(K,M,a);

    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

    [Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
    S013_Divided_G0=0.1;
    S023_Divided_G0=0;
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
    log_G1DG0(np,1)=-3+del*np;
    GeDG0(np,1)=delta_cof(1,1);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

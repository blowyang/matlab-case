clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
%在确定的模量比例下，和确定的界面相厚度的前提下，改变纤维的体积分数VF2,研究
%等效切边模量的变化规律

R2=1;

G0=100;
G2=1;

R1=R2+0.2;



N=30;
M=25;
%k层界面相
K=1;
%KP个配点
KP=50;


NP=40;
%纤维体积分数增量
del=0.01;
vf2=zeros(NP,1);
GeDG0=zeros(NP,1);
for np=1:NP
    VF2=del*np;
    a=(pi/VF2)^0.5*R2;
    lambda1=R1/a;
    
    [E_Matrix,N_Matrix]=Get_E_N_Matrix(N,K,R1,R2,G0,G2);
    B_Matrix=Get_B_Matrix(N,M,K,R1,R2,a,G0,G2);
    C_Matrix=Get_C_Matrix(KP,N,a,R1);
    D_Matrix=Get_D_Matrix(KP,M,a);

    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

    [Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(KP,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
    S013_Divided_G0=0.1;
    S023_Divided_G0=0;
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
    vf2(np,1)=del*np;
    GeDG0(np,1)=delta_cof(1,1);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

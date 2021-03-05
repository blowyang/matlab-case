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

x1=-1;
x2=1;
G0=1000;
G1=10;
G2=10;
a=4;

N=30;
M=25;
K=40;


NP=16;
%纤维体积分数增量
del=0.05;
R0=0.1;
GeDG0=zeros(NP,3);
for np=1:NP
    R=R0+del*np;
    GeDG0(np,1)=2*pi*R^2/a^2;
    
    R1=R;
    R2=R;
    
    [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,R1,R2,x1,x2,G0,G1,G2);
    B_Matrix=Get_B_Matrix(N,M,x1,R1,x2,R2,a);
    C_Matrix=Get_C_Matrix(K,N,a,R1,R2,x1,x2);
    D_Matrix=Get_D_Matrix(K,M,a);
    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
    [Xcd1,Xcd2,Xe1,Xe2]=Get_Xc_Xe_SubElement(K,H_Matrix,B_Matrix,F_Matrix);
    S013_Divided_G0=0.1;
    S023_Divided_G0=0;
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,R2,x1,x2,a,Xcd1,Xcd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
    
    GeDG0(np,2)=delta_cof(1,1);
    GeDG0(np,3)=delta_cof(2,2);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

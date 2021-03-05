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

x1=-10;
x2=10;
G0=1000;
G1=10;
G2=10;
a=40;

N=30;
M=25;
K=40;

R1=1;
R2=1;

[~,F_Matrix]=Get_E_F_Matrix(N,R1,R2,x1,x2,G0,G1,G2);
B_Matrix=Get_B_Matrix(N,M,x1,R1,x2,R2,a);
C_Matrix=Get_C_Matrix(K,N,a,R1,R2,x1,x2);
D_Matrix=Get_D_Matrix(K,M,a);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
[Xcd1,Xcd2,Xe1,Xe2]=Get_Xc_Xe_SubElement(K,H_Matrix,B_Matrix,F_Matrix);

S013_Divided_G0=0;
S023_Divided_G0=1/3^0.5;

[delta1_Divided_a,delta2_Divided_a,~]=Get_Delta(R1,R2,x1,x2,a,Xcd1,Xcd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);

delta1=delta1_Divided_a*a;
delta2=delta2_Divided_a*a;
Xcd=delta1*Xcd1+delta2*Xcd2;
Xe=delta1*Xe1+delta2*Xe2;

COU=61;
del_2=2*pi/(COU+1);
mises_L1=zeros();
for cou=1:COU
    zM_np=x2+R2*cos(del_2*cou)+R2*sin(del_2*cou)*1i;
    MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_np); 
    mises_L1(cou,1)=del_2*cou;
    mises_L1(cou,2)=MisesMatrix;
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
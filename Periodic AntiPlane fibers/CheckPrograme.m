clear;
clc;

% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.

N=5;
M=5;
K=30;

R1=0.6;
R2=0.6;
x1=-1;
x2=1;
G0=10;
G1=100;
G2=100;
a=4;

[E_Matrix,F_Matrix]=Get_E_F_Matrix(N,R1,R2,x1,x2,G0,G1,G2);
B_Matrix=Get_B_Matrix(N,M,x1,R1,x2,R2,a);
C_Matrix=Get_C_Matrix(K,N,a,R1,R2,x1,x2);
D_Matrix=Get_D_Matrix(K,M,a);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
[Xcd1,Xcd2,Xe1,Xe2]=Get_Xc_Xe_SubElement(K,H_Matrix,B_Matrix,F_Matrix);
S013_Divided_G0=0;
S023_Divided_G0=0.1;
[delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,R2,x1,x2,a,Xcd1,Xcd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);

delta1=delta1_Divided_a*a;
delta2=delta2_Divided_a*a;
Xcd=delta1*Xcd1+delta2*Xcd2;
Xe=delta1*Xe1+delta2*Xe2;

NP=100;
del=a/(NP+1);
w_AB=zeros(NP,1);
w_CD=zeros(NP,1);
w_AB_Sub_CD=zeros(NP,1);
w_DB=zeros(NP,1);
w_CA=zeros(NP,1);
w_DB_Sub_CA=zeros(NP,1);
S23_AB=zeros(NP,1);
S13_DB=zeros(NP,1);
S13_CA=zeros(NP,1);
S23_CD=zeros(NP,1);
S13_DB_Sub_CA=zeros(NP,1);
S23_AB_Sub_CD=zeros(NP,1);
for np=1:NP
    zk_AB=del*np-a/2+a*1i/2;
    zk_DB=del*np*1i-a*1i/2+a/2;
    zk_CD=conj(zk_AB);
    zk_CA=conj(-zk_DB);
    
    w_AB(np,1)=Get_W0_Dispalcement(a,R1,R2,x1,x2,Xcd,Xe,zk_AB); 
    w_CD(np,1)=Get_W0_Dispalcement(a,R1,R2,x1,x2,Xcd,Xe,zk_CD);
    w_CA(np,1)=Get_W0_Dispalcement(a,R1,R2,x1,x2,Xcd,Xe,zk_CA); 
    w_DB(np,1)=Get_W0_Dispalcement(a,R1,R2,x1,x2,Xcd,Xe,zk_DB);
    
    w_AB_Sub_CD(np,1)=(w_AB(np,1)-w_CD(np,1))/(delta2);
    w_DB_Sub_CA(np,1)=(w_DB(np,1)-w_CA(np,1))/(delta1);
    
    S23_AB(np,1)=Get_S023_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zk_AB);
    S23_CD(np,1)=Get_S023_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zk_CD);
    S13_DB(np,1)=Get_S013_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zk_DB);
    S13_CA(np,1)=Get_S013_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zk_CA);
 
    S13_DB_Sub_CA(np,1)=S13_DB(np,1)-S13_CA(np,1);
    S23_AB_Sub_CD(np,1)=S23_AB(np,1)-S23_CD(np,1);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

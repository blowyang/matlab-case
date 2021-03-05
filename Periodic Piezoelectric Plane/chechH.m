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
K=40;
Gm=10;
Gf=100;
kf=2;
km=2;
R=1;

a=4;
e31=-6.5;
%E3infinity=0.2;
E3infinity=0;
[E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km);
B_Matrix=Get_B_Matrix(N,M,km,R,a);
ELoad=Get_ELoad(N,e31,E3infinity,R);
[C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,a,R,km,Gm);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
[Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,C_Matrix,F_Matrix,ELoad);
S022=0;
S011=0;
S012=1;
% [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
%     Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
[delta1,delta2,delta3,delta4,delta_cof]=Get_Delta(a,R,Xce0,Xce1,Xce2,Xce3,Xce4,...
    Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
Xce=Xce0+delta1*Xce1+delta2*Xce2+delta3*Xce3+delta4*Xce4;
Xdf=Xdf0+delta1*Xdf1+delta2*Xdf2+delta3*Xdf3+delta4*Xdf4;


NP=100;
del=a/(NP+1);

u2_AB=zeros(NP,1);
u2_CD=zeros(NP,1);
u2_AB_Sub_CD=zeros(NP,1);
u1_AB=zeros(NP,1);
u1_CD=zeros(NP,1);
u1_AB_Sub_CD=zeros(NP,1);
% w_DB=zeros(NP,1);
% w_CA=zeros(NP,1);
% w_DB_Sub_CA=zeros(NP,1);
% S23_AB=zeros(NP,1);
% S13_DB=zeros(NP,1);
% S13_CA=zeros(NP,1);
% S23_CD=zeros(NP,1);
% S13_DB_Sub_CA=zeros(NP,1);
% S23_AB_Sub_CD=zeros(NP,1);
for np=1:NP
    z_AB=del*np-a/2+a*1i/2;
    z_DB=del*np*1i-a*1i/2+a/2;
    z_CD=conj(z_AB);
    z_CA=conj(-z_DB);
    
    [u1,u2]=Get_Dispalcement(a,R,Gm,km,Xce,Xdf,z_AB);
    u1_AB(np,1)=u1;
    u2_AB(np,1)=u2;
    [u1,u2]=Get_Dispalcement(a,R,Gm,km,Xce,Xdf,z_CD);
    u1_CD(np,1)=u1;
    u2_CD(np,1)=u2;
    u2_AB_Sub_CD(np,1)=(u2_AB(np,1)-u2_CD(np,1))/(delta2);
    u1_AB_Sub_CD(np,1)=(u1_AB(np,1)-u1_CD(np,1))/(delta1);
%     w_CA(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_CA); 
%     w_DB(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_DB);
%     
%     w_AB_Sub_CD(np,1)=(w_AB(np,1)-w_CD(np,1))/(delta2);
%     w_DB_Sub_CA(np,1)=(w_DB(np,1)-w_CA(np,1))/(delta1);
%     
%     S23_AB(np,1)=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,z_AB);
%     S23_CD(np,1)=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,z_CD);
%     S13_DB(np,1)=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,z_DB);
%     S13_CA(np,1)=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,z_CA);
%  
%     S13_DB_Sub_CA(np,1)=S13_DB(np,1)-S13_CA(np,1);
%     S23_AB_Sub_CD(np,1)=S23_AB(np,1)-S23_CD(np,1);
end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.
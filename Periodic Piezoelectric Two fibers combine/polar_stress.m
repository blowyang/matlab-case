clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
N=30;
M=20;
K=100;
Gm=10;
Gf1=100;
kf1=1.8;
Gf2=1;
kf2=1.8;
km=1.8;
R1=0.8;
R2=0.8;
x1=-0.9;
x2=0.9;
e31_1=-6.5;
E3infinity_1=0;
e31_2=-6.5;
E3infinity_2=0.03;
VF=0.3;
a=sqrt(pi*(R1^2+R2^2)/VF);
%e31=-6.5;
%E3infinity=0.2;
%E3infinity=0;
A_Matrix=Get_A_Matrix(N,R1,R2,x1,x2,Gm,km,Gf1,kf1,Gf2,kf2);
B_Matrix=Get_B_Matrix(N,M,R1,R2,x1,x2,a,km);
[E_Matrix,F_Matrix]=Get_E_F_Matrix(N,A_Matrix);
[C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,R1,R2,x1,x2,a,km,Gm);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
ELoad=Get_ELoad(N,e31_1,E3infinity_1,R1,e31_2,E3infinity_2,R2);
%[Xce1,Xce2,Xce3,Xce4,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix);
[Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix,C_Matrix,ELoad);
S022=0;
S011=1;
S012=0;
% [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce1,Xce2,Xce3,Xce4,...
%     Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
[delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce0,Xce1,Xce2,Xce3,Xce4,...
    Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
% 
Xce=Xce0+delta1*Xce1+delta2*Xce2+delta3*Xce3+delta4*Xce4;
Xdf=Xdf0+delta1*Xdf1+delta2*Xdf2+delta3*Xdf3+delta4*Xdf4;


NP=100;
del=2*pi/(NP+1);
stress=zeros(NP,3);
for np=1:NP
    theta=np*del;
    z=x2+R2*exp(1i*theta);
    
    [Sm22,Sm11,Sm12]=Get_Matrix_Stress(R1,R2,x1,x2,a,Xce,Xdf,z);
    
    Smrr=(Sm22+Sm11)/2+(Sm11-Sm22)/2*cos(2*theta)+Sm12*sin(2*theta);
    Smss=(Sm22+Sm11)/2-(Sm11-Sm22)/2*cos(2*theta)-Sm12*sin(2*theta);
    stress(np,1)=theta;
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
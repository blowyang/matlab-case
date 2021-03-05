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
K=80;
Gm=1;
Gf1=Gm;
kf1=1.8;
Gf2=10;
kf2=1.8;
km=1.8;
a=4;
R1=0.1;
x1=-1.75;
x2=0;
e31_1=-6.5;
E3infinity_1=0;
e31_2=-6.5;
E3infinity_2=0;

S022=1;
S011=0;
S012=0;
% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
Em=2.6*Gm;

NP=30;
del=0.012;
Stif=zeros(NP,8);
for np=1:NP
    VF=del*np;
    
    R2=sqrt(VF/pi)*a;
    
    A_Matrix=Get_A_Matrix(N,R1,R2,x1,x2,Gm,km,Gf1,kf1,Gf2,kf2);
    B_Matrix=Get_B_Matrix(N,M,R1,R2,x1,x2,a,km);
    [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,A_Matrix);
    [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,R1,R2,x1,x2,a,km,Gm);
    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
    ELoad=Get_ELoad(N,e31_1,E3infinity_1,R1,e31_2,E3infinity_2,R2);
    [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix,C_Matrix,ELoad);
    
    [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce0,Xce1,Xce2,Xce3,Xce4,...
        Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
    
    Stif(np,1)=VF;
    
    %S11,S22
    tem=delta_cof(1,2);
    Stif(np,2)=tem*a/Em;
    %S12
    tem=delta_cof(1,3);
    Stif(np,3)=tem*a/Em;
    %S33
    tem=delta_cof(2,1);
    Stif(np,4)=tem*a/Gm;
    %S12,S21,S13,S23
    tem=delta_cof(2,2);
    Stif(np,5)=tem*a/Gm;
    tem=delta_cof(2,3);
    Stif(np,6)=tem*a/Gm;
    tem=delta_cof(4,2);
    Stif(np,7)=tem*a/Gm;
    tem=delta_cof(4,3);
    Stif(np,8)=tem*a/Gm;
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
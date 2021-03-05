clear;
clc;
% 启动并行计算
core_number=10;            %想要调用的处理器个数
parpool('local',core_number);
% % % % % % 启动后有如下提示：
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
N=25;
M=18;
K=80;

a=4;
C11E=12.6;
C12E=5.5;
Ef=C11E-C12E/C11E;
muf=C12E/C11E;

Gf1=Ef/2/(1+muf);
kf1=3-4*muf;

Gm=5*Gf1/10;
km=1.8;

Gf2=Gf1;
kf2=kf1;

x1=-(2^0.5-1)/2*a-(2^0.5-1)/2*a*1i;
x2=(2^0.5-1)/2*a+(2^0.5-1)/2*a*1i;
% x1=-1;
% x2=1;

e31_1=-6.5;
E3infinity_1=1;
e31_2=-6.5;
E3infinity_2=1;

S022=0;
S011=0;
S012=0;

% VF=0.4;
% 
% R1=sqrt(VF/2/pi)*a;
% R2=R1;
% 
% [A_Matrix,B_Matrix]=Get_AB_Matrix(N,M,R1,R2,x1,x2,a,Gm,km,Gf1,kf1,Gf2,kf2);
% [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,A_Matrix);
% [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,R1,R2,x1,x2,a,km,Gm);
% H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
% ELoad=Get_ELoad(N,e31_1,E3infinity_1,R1,e31_2,E3infinity_2,R2);
% [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix,C_Matrix,ELoad);
% 
% [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce0,Xce1,Xce2,Xce3,Xce4,...
%     Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);


% delta1=1;
% delta2=1;
% delta3=1;
% delta4=1;
Em=2.6*Gm;

NP=30;
del=0.012;
Stif=zeros(NP,3);
for np=1:NP
    VF=del*np;
    
    R1=sqrt(VF/2/pi)*a;
    R2=R1;
    [A_Matrix,B_Matrix]=Get_AB_Matrix(N,M,R1,R2,x1,x2,a,Gm,km,Gf1,kf1,Gf2,kf2);
    [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,A_Matrix);
    [C_Matrix,D_Matrix]=Get_CD_Matrix(K,N,M,R1,R2,x1,x2,a,km,Gm);
    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,F_Matrix,B_Matrix);
    ELoad=Get_ELoad(N,e31_1,E3infinity_1,R1,e31_2,E3infinity_2,R2);
    [Xce0,Xce1,Xce2,Xce3,Xce4,Xdf0,Xdf1,Xdf2,Xdf3,Xdf4]=Get_Xce_Xdf_SubElement(K,H_Matrix,B_Matrix,F_Matrix,C_Matrix,ELoad);
   
    [delta1,delta2,delta3,delta4,delta_cof]=Solve_Delta(R1,R2,x1,x2,a,Xce0,Xce1,Xce2,Xce3,Xce4,...
        Xdf0,Xdf1,Xdf2,Xdf3,Xdf4,S022,S011,S012);
    
    Stif(np,1)=VF;
    Stif(np,2)=delta3/a;
    Stif(np,3)=delta2/a;
    
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
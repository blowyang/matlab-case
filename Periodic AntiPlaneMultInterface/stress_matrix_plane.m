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
x1=0;

R1=R2+0.2;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=100;
G2=1;

N=30;
M=25;
%K 层数
%K=25;
K=25;
%KP 配点数
KP=20;


mises_matrix=zeros();
mises_L1=zeros();
%纤维体积分数增量
S023_Divided_G0=1/6^0.5;
S013_Divided_G0=1/6^0.5;
[E_Matrix,N_Matrix]=Get_E_N_Matrix(N,K,R1,R2,G0,G2);
B_Matrix=Get_B_Matrix(N,M,K,R1,R2,a,G0,G2);
C_Matrix=Get_C_Matrix(KP,N,a,R1);
D_Matrix=Get_D_Matrix(KP,M,a);
H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

[Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(KP,H_Matrix,B_Matrix,E_Matrix,N_Matrix);

[delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
delta1=delta1_Divided_a*a;
delta2=delta2_Divided_a*a;
Xa=delta1*Xa2+delta2*Xa1;
Xd=delta1*Xd2+delta2*Xd1;
Xe=delta1*Xe2+delta2*Xe1;

% n
n=1;

NR=10;
dr_1=R1/100;

COU=60;
del_2=2*pi/(COU+1);
for nr=1:NR
    R1_Show=R1+dr_1*(nr-1);
    for cou=1:COU+1
        theta=del_2*cou;
        zM_circle_1=x1+R1_Show*exp(theta*1i);
        MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_circle_1);
        mises_matrix(n,1)=real(zM_circle_1);
        mises_matrix(n,2)=imag(zM_circle_1);
        mises_matrix(n,3)=MisesMatrix;
        n=n+1; 
    end    
end

NP=60;
del_1=a/(NP+1);

for np_x=1:NP  
    zx=-a/2+del_1*np_x;
    for np_y=1:NP
        zy=-a/2+del_1*np_y;
        if (abs(zx+zy*1i-x1)>=R1+dr_1*(NR-1) )
            zM_np=zx+zy*1i;
            MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_np);
            mises_matrix(n,1)=real(zM_np);
            mises_matrix(n,2)=imag(zM_np);
            mises_matrix(n,3)=MisesMatrix;
            n=n+1;
        end       
    end
end

NP_B=80;
del_b=a/(NP_B+1);
for np=1:NP_B
    zM_left=-a/2-a/2*1i+del_b*np*1i;
    MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_left);
    mises_matrix(n,1)=real(zM_left);
    mises_matrix(n,2)=imag(zM_left);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_right=a/2-a/2*1i+del_b*np*1i;
    MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_right);
    mises_matrix(n,1)=real(zM_right);
    mises_matrix(n,2)=imag(zM_right);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_top=a/2*1i-a/2+del_b*np;
    MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_top);
    mises_matrix(n,1)=real(zM_top);
    mises_matrix(n,2)=imag(zM_top);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_bottom=-a/2*1i-a/2+del_b*np;
    MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_bottom);
    mises_matrix(n,1)=real(zM_bottom);
    mises_matrix(n,2)=imag(zM_bottom);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;    
end

% four conner points
zM_conner_1=-a/2-a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_conner_1);
mises_matrix(n,1)=real(zM_conner_1);
mises_matrix(n,2)=imag(zM_conner_1);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_2=a/2-a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_conner_2);
mises_matrix(n,1)=real(zM_conner_2);
mises_matrix(n,2)=imag(zM_conner_2);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_3=-a/2+a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_conner_3);
mises_matrix(n,1)=real(zM_conner_3);
mises_matrix(n,2)=imag(zM_conner_3);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_4=a/2+a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_conner_4);
mises_matrix(n,1)=real(zM_conner_4);
mises_matrix(n,2)=imag(zM_conner_4);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
x = mises_matrix(:,1); 
y = mises_matrix(:,2); 
z = mises_matrix(:,3);
%title(['G_{m}:G_{f}=',num2str(Gm),':',num2str(Gf),'; e_{31}E_{0}^{\infty}:\sigma_{11}^{\infty}=',num2str(de/S011)]);%添加标题
shp=alphaShape(x,y,'HoleThreshold',0.5);
tri = alphaTriangulation(shp);

patch('faces',tri,'vertices',[x,y],... 
    'facevertexcdata',z,'edgecolor','none','facecolor','interp'); 


axis image
colorbar EastOutside

% for cou=1:COU
%     zM_np=R1*cos(del_2*cou)+R1*sin(del_2*cou)*1i;
%     MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_np); 
%     mises_L1(cou,1)=del_2*cou;
%     mises_L1(cou,2)=MisesMatrix;
% end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

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


G0=10;
G1=1;
G2=1;
a=4;
x1=-(0.207*a+0.207*a*1i);
x2=(0.207*a+0.207*a*1i);
% x1=-1;
% x2=1;
N=30;
M=25;
K=40;

R1=0.8;
R2=0.8;

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

mises_matrix=zeros();
mises_L1=zeros();
%纤维体积分数增量
% n
n=1;

NR=10;
dr_1=R1/100;
dr_2=R2/100;

COU=60;
del_2=2*pi/(COU+1);
for nr=1:NR
    R1_Show=R1+dr_1*(nr-1);
    R2_Show=R2+dr_2*(nr-1);
    for cou=1:COU+1
        theta=del_2*cou;
        zM_circle_1=x1+R1_Show*exp(theta*1i);
        MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_circle_1);
        mises_matrix(n,1)=real(zM_circle_1);
        mises_matrix(n,2)=imag(zM_circle_1);
        mises_matrix(n,3)=MisesMatrix;
        n=n+1; 
        zM_circle_2=x2+R2_Show*exp(theta*1i);
        MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_circle_2);
        mises_matrix(n,1)=real(zM_circle_2);
        mises_matrix(n,2)=imag(zM_circle_2);
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
        if (abs(zx+zy*1i-x1)>=R1+dr_1*(NR-1) ) && (abs(zx+zy*1i-x2)>=R2+dr_1*(NR-1) )
            zM_np=zx+zy*1i;
            MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_np);
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
    MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_left);
    mises_matrix(n,1)=real(zM_left);
    mises_matrix(n,2)=imag(zM_left);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_right=a/2-a/2*1i+del_b*np*1i;
    MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_right);
    mises_matrix(n,1)=real(zM_right);
    mises_matrix(n,2)=imag(zM_right);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_top=a/2*1i-a/2+del_b*np;
    MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_top);
    mises_matrix(n,1)=real(zM_top);
    mises_matrix(n,2)=imag(zM_top);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;
    zM_bottom=-a/2*1i-a/2+del_b*np;
    MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_bottom);
    mises_matrix(n,1)=real(zM_bottom);
    mises_matrix(n,2)=imag(zM_bottom);
    mises_matrix(n,3)=MisesMatrix;
    n=n+1;    
end

% four conner points
zM_conner_1=-a/2-a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_conner_1);
mises_matrix(n,1)=real(zM_conner_1);
mises_matrix(n,2)=imag(zM_conner_1);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_2=a/2-a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_conner_2);
mises_matrix(n,1)=real(zM_conner_2);
mises_matrix(n,2)=imag(zM_conner_2);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_3=-a/2+a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_conner_3);
mises_matrix(n,1)=real(zM_conner_3);
mises_matrix(n,2)=imag(zM_conner_3);
mises_matrix(n,3)=MisesMatrix;
n=n+1;
zM_conner_4=a/2+a/2*1i;
MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_conner_4);
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

% NPP=120;
% del_a=2*pi/NPP;
% 
% for npp=1:NPP
%     theta=del_a*(npp-1);
%     zM_npp=x1+R1*cos(theta)+R1*sin(theta)*1i;
%     MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_npp);
%     mises_L1(npp,1)=theta;
%     mises_L1(npp,2)=MisesMatrix;
% end

%% 关闭并行计算
delete(gcp('nocreate'));
% % % % % 关闭后有如下提示：
% Parallel pool using the 'local' profile is shutting down.

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
G0=1;
G1=10;
G2=10;
a=4;

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
NP=70;
del_1=0.75*a/(NP+1);
COU=60;
del_2=2*pi/(COU+1);

for np=1:NP+2
    R=del_1*(np);
    
    for cou=1:COU+1
        zx=R*cos(del_2*cou);
        zy=R*sin(del_2*cou);
        if(abs(zx)<=a/2) && (abs(zy)<=a/2) && (abs(zx+zy*1i-1)>=0.8 ) && (abs(zx+zy*1i+1)>=0.8 )
            zM_np=zx+zy*1i;
            MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,zM_np);
            mises_matrix(n,1)=R*cos(del_2*cou);
            mises_matrix(n,2)=R*sin(del_2*cou);
            mises_matrix(n,3)=MisesMatrix;
            n=n+1;
        end
        
    end
end

x = mises_matrix(:,1); 
y = mises_matrix(:,2); 
z = mises_matrix(:,3);

shp=alphaShape(x,y,'HoleThreshold',1.5);
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

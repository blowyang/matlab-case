clear;
clc;

R2=1;
R1=R2+0.2;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=100;
G1=10;
G2=1;
gama1=G0/G1;
gama2=G1/G2;
lambda=R2/R1;
lambda1=R1/a;

N=30;
M=25;
%K 配点数
K=20;

NP=15;
%t13的增量
del=1/3^0.5/(NP+1);

mises_max=zeros(2*NP,2);

for np=1:NP
    
    S013_Divided_G0=del*np;
    S023_Divided_G0=(1/3-S013_Divided_G0^2)^0.5;
    
    [E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
    B_Matrix=Get_B_Matrix(N,M,gama1,lambda1);

    C_Matrix=Get_C_Matrix(K,N,a,R1);
    D_Matrix=Get_D_Matrix(K,M,a);

    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

    [Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
    
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,-S023_Divided_G0);
    delta1=delta1_Divided_a*a;
    delta2=delta2_Divided_a*a;
    Xa=delta1*Xa2+delta2*Xa1;
    Xd=delta1*Xd2+delta2*Xd1;
    Xe=delta1*Xe2+delta2*Xe1;
    COU=100;
    DELT=2*pi/(COU+1);
    last_max=0;
    for cou=1:COU+1
        zM_np=R1*cos(DELT*cou)+R1*sin(DELT*cou)*1i;
        MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_np);
        if (MisesMatrix>last_max)
            last_max=MisesMatrix;
        end          
    end
    mises_max(2*np-1,1)=-S013_Divided_G0/S023_Divided_G0;
    mises_max(2*np-1,2)=last_max;
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
    delta1=delta1_Divided_a*a;
    delta2=delta2_Divided_a*a;
    Xa=delta1*Xa2+delta2*Xa1;
    Xd=delta1*Xd2+delta2*Xd1;
    Xe=delta1*Xe2+delta2*Xe1;
    COU=100;
    DELT=2*pi/(COU+1);
    last_max=0;
    for cou=1:COU+1
        zM_np=R1*cos(DELT*cou)+R1*sin(DELT*cou)*1i;
        MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,zM_np);
        if (MisesMatrix>last_max)
            
            last_max=MisesMatrix;
        end          
    end
    mises_max(2*np,1)=S013_Divided_G0/S023_Divided_G0;
    mises_max(2*np,2)=last_max;
end

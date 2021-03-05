clear;
clc;
% �������м���
core_number=10;            %��Ҫ���õĴ���������
parpool('local',core_number);
% % % % % % ��������������ʾ��
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 2 workers.
%��ȷ����ģ�������£���ȷ���Ľ������ȵ�ǰ���£��ı���ά���������VF2,�о�
%��Ч�б�ģ���ı仯����

R2=1;
R1=R2+0.2;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=100;
G2=1;

N=30;
M=25;
%K ����
%K=25;
K=1;
%KP �����
KP=20;

[E_Matrix,N_Matrix]=Get_E_N_Matrix(N,K,R1,R2,G0,G2);

B_Matrix=Get_B_Matrix(N,M,K,R1,R2,a,G0,G2);
C_Matrix=Get_C_Matrix(KP,N,a,R1);
D_Matrix=Get_D_Matrix(KP,M,a);

H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

[Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(KP,H_Matrix,B_Matrix,E_Matrix,N_Matrix);

NP=15;
%��ά�����������
del=1/3^0.5/(NP+1);

mises_max=zeros(2*NP,2);

for np=1:NP
    
    S013_Divided_G0=del*np;
    S023_Divided_G0=(1/3-S013_Divided_G0^2)^0.5;
    
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

%% �رղ��м���
delete(gcp('nocreate'));
% % % % % �رպ���������ʾ��
% Parallel pool using the 'local' profile is shutting down.

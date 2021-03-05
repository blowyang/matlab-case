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

G0=100;
G1=10;
G2=1;

R1=R2+0.1;
gama1=G0/G1;
gama2=G1/G2;
lambda=R2/R1;
VF2=0.2;
a=(pi/VF2)^0.5*R2;


N=30;
M=25;
K=20;


NP=20;
%��ά�����������
del=0.01;
ds1=zeros(NP,1);
GeDG0=zeros(NP,1);
for np=1:NP
    R1=R2+del*np;
    lambda1=R1/a;
    [E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
    B_Matrix=Get_B_Matrix(N,M,gama1,lambda1);

    C_Matrix=Get_C_Matrix(K,N,a,R1);
    D_Matrix=Get_D_Matrix(K,M,a);

    H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

    [Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
    S013_Divided_G0=0.1;
    S023_Divided_G0=0;
    [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
    ds1(np,1)=del*np;
    GeDG0(np,1)=delta_cof(1,1);
end

%% �رղ��м���
delete(gcp('nocreate'));
% % % % % �رպ���������ʾ��
% Parallel pool using the 'local' profile is shutting down.

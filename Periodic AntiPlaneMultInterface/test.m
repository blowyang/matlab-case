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

G0=1;
G2=100;

R1=R2+0.2;



N=30;
M=25;
%k�������
K=20;
%KP�����
KP=50;

VF2=0.15;
a=(pi/VF2)^0.5*R2;
lambda1=R1/a;

[E_Matrix,N_Matrix]=Get_E_N_Matrix(N,K,R1,R2,G0,G2);
B_Matrix=Get_B_Matrix(N,M,K,R1,R2,a,G0,G2);
C_Matrix=Get_C_Matrix(KP,N,a,R1);
D_Matrix=Get_D_Matrix(KP,M,a);

H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

[Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(KP,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
S013_Divided_G0=1;
S023_Divided_G0=1;
[delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
%% �رղ��м���
delete(gcp('nocreate'));
% % % % % �رպ���������ʾ��
% Parallel pool using the 'local' profile is shutting down.
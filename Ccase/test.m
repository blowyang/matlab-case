clear;
clc;
R2=1;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=1250;
gama1=0.1;
G1=G0/gama1;
gama2=0.1;
G2=G1/gama2;
R1=R2+0.1;
lambda=R2/R1;
lambda1=R1/a ;

N=30;
K=20;
[E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
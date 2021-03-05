function [E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2)
A_Matrix=Get_A_Matrix(N,lambda,gama1,gama2);
A_Matrix=inv(A_Matrix);
E_Matrix=zeros(2*N,8*N);
F_Matrix=zeros(2*N,8*N);
G_Matrix=zeros(2*N,8*N);
N_Matrix=zeros(2*N,8*N);
for n=1:2*N
    E_Matrix(n,:)=A_Matrix(n,:);
    F_Matrix(n,:)=A_Matrix(n+2*N,:);
    G_Matrix(n,:)=A_Matrix(n+4*N,:);
    N_Matrix(n,:)=A_Matrix(n+6*N,:);
end
end
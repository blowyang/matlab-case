function [E_Matrix,N_Matrix]=Get_E_N_Matrix(N,K,R1,R2,G0,G2)
A_Matrix=Get_A_Matrix(N,K,R1,R2,G0,G2);
A_Matrix=inv(A_Matrix);
E_Matrix=zeros(2*N,(2+2*K)*2*N);
N_Matrix=zeros(2*N,(2+2*K)*2*N);
for n=1:2*N
    E_Matrix(n,:)=A_Matrix(n,:);
    N_Matrix(n,:)=A_Matrix(4*K*N+2*N+n,:);
end
end
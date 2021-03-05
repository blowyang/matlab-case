function [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,A_Matrix)
A_Matrix=inv(A_Matrix);
E_Matrix=zeros(8*N,16*N);
F_Matrix=zeros(8*N,16*N);

for n=1:8*N
    E_Matrix(n,:)=A_Matrix(n,:);
    F_Matrix(n,:)=A_Matrix(n+8*N,:);
end
end
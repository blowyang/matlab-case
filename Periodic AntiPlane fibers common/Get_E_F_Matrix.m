function [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,R1,R2,x1,x2,G0,G1,G2)
A_Matrix=Get_A_Matrix(N,R1,R2,x1,x2,G0,G1,G2);
A_Matrix=inv(A_Matrix);
E_Matrix=zeros(4*N,8*N);
F_Matrix=zeros(4*N,8*N);

for n=1:4*N
    E_Matrix(n,:)=A_Matrix(n,:);
    F_Matrix(n,:)=A_Matrix(n+4*N,:);
end
end
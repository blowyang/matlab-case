function [E_Matrix,F_Matrix]=Get_E_F_Matrix(N,Gm,Gf,kf,km)
A_Matrix=Get_A_Matrix(N,Gm,Gf,kf,km);
A_Matrix=inv(A_Matrix);
E_Matrix=zeros(4*N,8*N);
F_Matrix=zeros(4*N,8*N);

for n=1:4*N
    E_Matrix(n,:)=A_Matrix(n,:);
    F_Matrix(n,:)=A_Matrix(n+4*N,:);
end
end
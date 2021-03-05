function B_Matrix=Get_B_Matrix(N,lambda1,gama1)
%lambda1=R1/a
%gama1=G0/G1
B_Matrix=zeros(8*N,2*N);
for i=1:N
    B_Matrix(i+4*N,2*i-1)=-lambda1^i;
    B_Matrix(i+5*N,2*i)=-lambda1^i;
    B_Matrix(i+6*N,2*i-1)=-gama1*lambda1^i;
    B_Matrix(i+7*N,2*i)=-gama1*lambda1^i;
end
end
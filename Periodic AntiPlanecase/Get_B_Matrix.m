function B_Matrix=Get_B_Matrix(N,M,gama1,lambda1)
%lambda1=R1/a 
%gama1=G0/G1
B_Matrix=zeros(8*N,2*M);
for i=1:N
    for j=1:M
        B_Matrix(i+4*N,2*j-1)=-Get_P_Coef(lambda1,j,i);
        B_Matrix(i+5*N,2*j)=-Get_P_Coef(lambda1,j,i);
        B_Matrix(i+6*N,2*j-1)=-gama1*Get_P_Coef(lambda1,j,i);
        B_Matrix(i+7*N,2*j)=-gama1*Get_P_Coef(lambda1,j,i);
    end    
end
end
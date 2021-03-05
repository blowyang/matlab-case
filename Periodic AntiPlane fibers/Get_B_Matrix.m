function B_Matrix=Get_B_Matrix(N,M,x1,R1,x2,R2,a)
%lambda1=R1/a 
%gama1=G0/G1
B_Matrix=zeros(8*N,2*M);
for j=1:N
    for k=1:M
        B_Matrix(j,2*k-1)=-Get_P_Coef(k,j,x1,R1,a);
        B_Matrix(j+N,2*k)=-Get_P_Coef(k,j,x1,R1,a);
        B_Matrix(j+2*N,2*k-1)=-Get_P_Coef(k,j,x1,R1,a);
        B_Matrix(j+3*N,2*k)=-Get_P_Coef(k,j,x1,R1,a);
        B_Matrix(j+4*N,2*k-1)=-Get_P_Coef(k,j,x2,R2,a);
        B_Matrix(j+5*N,2*k)=-Get_P_Coef(k,j,x2,R2,a);
        B_Matrix(j+6*N,2*k-1)=-Get_P_Coef(k,j,x2,R2,a);
        B_Matrix(j+7*N,2*k)=-Get_P_Coef(k,j,x2,R2,a);
    end    
end
end
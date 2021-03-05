function B_Matrix=Get_B_Matrix(N,M,K,R1,R2,a,G0,G2)
%lambda1=R1/a 
%gama1=G0/G1
h=(R1-R2)/K;
B_Matrix=zeros((2+2*K)*2*N,2*M);
lambda1=R1/a;
Rk=R1-0.5*h;
gama1=G0/Get_G1(G0,G2,R1,R2,Rk);
for i=1:N
    for j=1:M
        B_Matrix(i+4*K*N,2*j-1)=-Get_P_Coef(lambda1,j,i);
        B_Matrix(i+4*K*N+N,2*j)=-Get_P_Coef(lambda1,j,i);
        B_Matrix(i+4*K*N+2*N,2*j-1)=-gama1*Get_P_Coef(lambda1,j,i);
        B_Matrix(i+4*K*N+3*N,2*j)=-gama1*Get_P_Coef(lambda1,j,i);
    end    
end
end
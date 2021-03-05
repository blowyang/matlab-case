function [B_Cof_Pos_j,B_Cof_Neg_j]=Get_B_Matrix_Cof_Three(N,R2,x2,a,j)

S0=zeros(N,2*N);
for n=1:N
    S0(n,2*n-1)=1;
    S0(n,2*n)=1i;
end

L2_Pos_j=zeros(1,N);
for n=1:N
   L2_Pos_j(1,n)=Cnm(n,j)*(x2/a)^(n-j)*(R2/a)^j;
end
L2_Neg_j=zeros(1,N);

zL4_Conj_Pos_j=zeros(1,N);
if j==1
    for n=1:N
        zL4_Conj_Pos_j(1,n)=R2*n/a*Cnm(n-1,0)*(x2/a)^(n-1);
    end    
end
zL4_Conj_Neg_One_j=zeros(1,N);
for n=1:N
    zL4_Conj_Neg_One_j(1,n)=x2*n/a*Cnm(n-1,j)*(x2/a)^(n-1-j)*(R2/a)^j;
end
zL4_Conj_Neg_Two_j=zeros(1,N);
for n=1:N
    zL4_Conj_Neg_Two_j(1,n)=R2*n/a*Cnm(n-1,j+1)*(x2/a)^(n-1-j-1)*(R2/a)^(j+1);
end
zL4_Conj_Neg_j=zL4_Conj_Neg_One_j+zL4_Conj_Neg_Two_j;

L2_Conj_Pos_j=zeros(1,N);
L2_Conj_Neg_j=zeros(1,N);
for n=1:N
   L2_Conj_Neg_j(1,n)=Cnm(n,j)*(x2/a)^(n-j)*(R2/a)^j;
end

B_Cof_Pos_j=[-L2_Pos_j*S0-zL4_Conj_Pos_j*conj(S0),-L2_Conj_Pos_j*conj(S0)];
B_Cof_Neg_j=[-L2_Neg_j*S0-zL4_Conj_Neg_j*conj(S0),-L2_Conj_Neg_j*conj(S0)];

end
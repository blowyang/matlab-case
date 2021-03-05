function [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_Three(N,R1,R2,x1,x2,j)

zero_list=zeros(1,2*N);

S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end

L72_Pos_j=zeros(1,N);
L72_Pos_j(1,j)=1;
L72_Neg_j=zeros(1,N);

zL82_Conj_Pos_j=zeros(1,N);
if j==1
    zL82_Conj_Pos_j(1,1)=1;   
end
zL82_Conj_Neg_j=zeros(1,N);
if(j+2<=N)
    zL82_Conj_Neg_j(1,j+2)=(j+2);%modifyed
end

L72_Conj_Pos_j=zeros(1,N);
L72_Conj_Neg_j=zeros(1,N);
L72_Conj_Neg_j(1,j)=1;

L11_Pos_j=zeros(1,N);
for n=1:N
    L11_Pos_j(1,n)=Cji(n,j)*(-1)^j*(R1/(x2-x1))^n*(R2/(x2-x1))^j;
end
L11_Neg_j=zeros(1,N);

zL31_Conj_Pos_j=zeros(1,N);
%added
if j==1
    for n=1:N
        zL31_Conj_Pos_j(1,n)=-R2*n/R1*(R1/(x2-x1))^(n+1);
    end
end
%added
zL31_Conj_Neg_One_j=zeros(1,N);
for n=1:N
    zL31_Conj_Neg_One_j(1,n)=-x2*n/R1*Cji(n+1,j)*(-1)^j*(R1/(x2-x1))^(n+1)*(R2/(x2-x1))^j;%modifyed
end
zL31_Conj_Neg_Two_j=zeros(1,N);
for n=1:N
    zL31_Conj_Neg_Two_j(1,n)=-R2*n/R1*Cji(n+1,j)*(-1)^(j+1)*(R1/(x2-x1))^(n+1)*(R2/(x2-x1))^(j+1);%modifyed
end
zL31_Conj_Neg_j=zL31_Conj_Neg_One_j+zL31_Conj_Neg_Two_j;

L12_Pos_j=zeros(1,N);
L12_Neg_j=zeros(1,N);
L12_Neg_j(1,j)=1;
%modifyed
zL32_Conj_Pos_One_j=zeros(1,N);
if j>=2
    zL32_Conj_Pos_One_j(1,j-1)=-x2*(j-1)/R2;   
end
zL32_Conj_Pos_Two_j=zeros(1,N);
if j>=3
    zL32_Conj_Pos_Two_j(1,j-2)=-(j-2);
end
zL32_Conj_Pos_j=zL32_Conj_Pos_One_j+zL32_Conj_Pos_Two_j;
zL32_Conj_Neg_j=zeros(1,N);
%modifyed
L12_Conj_Pos_j=zeros(1,N);
L12_Conj_Pos_j(1,j)=1;
L12_Conj_Neg_j=zeros(1,N);

L11_Conj_Pos_j=zeros(1,N);
L11_Conj_Neg_j=zeros(1,N);
for n=1:N
    L11_Pos_j(1,n)=Cji(n,j)*(-1)^j*(R1/(x2-x1))^n*(R2/(x2-x1))^j;
end
% disp("L71_Pos_j is");
% disp(L71_Pos_j*S1);
% disp("zL81_Conj_Pos_j is");
% disp(zL81_Conj_Pos_j*conj(S1));
A_Cof_Pos_j=[zero_list,zero_list,L72_Pos_j*S1+zL82_Conj_Pos_j*conj(S1),L72_Conj_Pos_j*conj(S1),...
    -L11_Pos_j*S1-zL31_Conj_Pos_j*conj(S1),-L12_Pos_j*S1-zL32_Conj_Pos_j*conj(S1),...
    -L11_Conj_Pos_j*conj(S1),-L12_Conj_Pos_j*conj(S1)];

A_Cof_Neg_j=[zero_list,zero_list,L72_Neg_j*S1+zL82_Conj_Neg_j*conj(S1),L72_Conj_Neg_j*conj(S1),...
    -L11_Neg_j*S1-zL31_Conj_Neg_j*conj(S1),-L12_Neg_j*S1-zL32_Conj_Neg_j*conj(S1),...
    -L11_Conj_Neg_j*conj(S1),-L12_Conj_Neg_j*conj(S1)];
end
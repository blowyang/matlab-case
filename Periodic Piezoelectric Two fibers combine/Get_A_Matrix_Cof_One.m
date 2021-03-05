function [A_Cof_Pos_j,A_Cof_Neg_j]=Get_A_Matrix_Cof_One(N,R1,R2,x1,x2,j)

zero_list=zeros(1,2*N);

S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end

L71_Pos_j=zeros(1,N);
L71_Pos_j(1,j)=1;
L71_Neg_j=zeros(1,N);

zL81_Conj_Pos_j=zeros(1,N);
if j==1
    zL81_Conj_Pos_j(1,1)=1;   
end
zL81_Conj_Neg_j=zeros(1,N);
if(j+2<=N)
    zL81_Conj_Neg_j(1,j+2)=(j+2);%modifyed
end

L71_Conj_Pos_j=zeros(1,N);
L71_Conj_Neg_j=zeros(1,N);
L71_Conj_Neg_j(1,j)=1;

L11_Pos_j=zeros(1,N);
L11_Neg_j=zeros(1,N);
L11_Neg_j(1,j)=1;

zL31_Conj_Pos_One_j=zeros(1,N);
if j>=2
    zL31_Conj_Pos_One_j(1,j-1)=-x1*(j-1)/R1;%modifyed
end
zL31_Conj_Pos_Two_j=zeros(1,N);
if j>=3
    zL31_Conj_Pos_Two_j(1,j-2)=-(j-2);%modifyed
end
zL31_Conj_Pos_j=zL31_Conj_Pos_One_j+zL31_Conj_Pos_Two_j;
zL31_Conj_Neg_j=zeros(1,N);

L12_Pos_j=zeros(1,N);
for n=1:N
    L12_Pos_j(1,n)=Cji(n,j)*(-1)^j*(R2/(x1-x2))^n*(R1/(x1-x2))^j;
end
L12_Neg_j=zeros(1,N);

zL32_Conj_Pos_j=zeros(1,N);
%added
if j==1
    for n=1:N
        zL32_Conj_Pos_j(1,n)=-R1*n/R2*(R2/(x1-x2))^(n+1);
    end
end
%added
zL32_Conj_Neg_One_j=zeros(1,N);
for n=1:N
    zL32_Conj_Neg_One_j(1,n)=-x1*n/R2*Cji(n+1,j)*(-1)^j*(R2/(x1-x2))^(n+1)*(R1/(x1-x2))^j;
end
zL32_Conj_Neg_Two_j=zeros(1,N);
for n=1:N
    zL32_Conj_Neg_Two_j(1,n)=-R1*n/R2*Cji(n+1,j)*(-1)^(j+1)*(R2/(x1-x2))^(n+1)*(R1/(x1-x2))^(j+1);
end
zL32_Conj_Neg_j=zL32_Conj_Neg_One_j+zL32_Conj_Neg_Two_j;

L11_Conj_Pos_j=zeros(1,N);
L11_Conj_Pos_j(1,j)=1;
L11_Conj_Neg_j=zeros(1,N);

L12_Conj_Pos_j=zeros(1,N);
L12_Conj_Neg_j=zeros(1,N);
for n=1:N
    L12_Pos_j(1,n)=Cji(n,j)*(-1)^j*(R2/(x1-x2))^n*(R1/(x1-x2))^j;
end
% disp("L71_Pos_j is");
% disp(L71_Pos_j*S1);
% disp("zL81_Conj_Pos_j is");
% disp(zL81_Conj_Pos_j*conj(S1));
A_Cof_Pos_j=[L71_Pos_j*S1+zL81_Conj_Pos_j*conj(S1),L71_Conj_Pos_j*conj(S1),zero_list,zero_list,...
    -L11_Pos_j*S1-zL31_Conj_Pos_j*conj(S1),-L12_Pos_j*S1-zL32_Conj_Pos_j*conj(S1),...
    -L11_Conj_Pos_j*conj(S1),-L12_Conj_Pos_j*conj(S1)];

A_Cof_Neg_j=[L71_Neg_j*S1+zL81_Conj_Neg_j*conj(S1),L71_Conj_Neg_j*conj(S1),zero_list,zero_list,...
    -L11_Neg_j*S1-zL31_Conj_Neg_j*conj(S1),-L12_Neg_j*S1-zL32_Conj_Neg_j*conj(S1),...
    -L11_Conj_Neg_j*conj(S1),-L12_Conj_Neg_j*conj(S1)];
end
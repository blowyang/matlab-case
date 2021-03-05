function [A_Cof,B_Cof]=AB_Matrix_Cof_One(N,M,R,a,z)
S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end
S0=zeros(M,2*M);
for m=1:M
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i;
end

[L1,L2,L3,L4,~,~,L7,L8]=Get_L_List(N,M,R,a,z);
A_Cof=[L7.'*S1+z*conj(L8.'*S1),conj(L7.'*S1),-(L1.'*S1+z*conj(L3.'*S1)),-conj(L1.'*S1)];
B_Cof=[-L2.'*S0-z*conj(L4.'*S0),-conj(L2.'*S0)];
end
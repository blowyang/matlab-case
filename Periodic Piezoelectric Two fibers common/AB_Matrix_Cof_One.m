function [A_Cof,B_Cof]=AB_Matrix_Cof_One(N,M,R1,R2,x1,x2,a,z)

zero_list=zeros(1,2*N);
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

[L11,L12,L2,L31,L32,L4,~,~,~,L71,~,L81,~]=Get_L_List(N,M,R1,R2,x1,x2,a,z);
A_Cof=[L71.'*S1+(z-x1)*conj(L81.'*S1),conj(L71.'*S1),zero_list,zero_list,...
    -(L11.'*S1+z*conj(L31.'*S1)),-(L12.'*S1+z*conj(L32.'*S1)),-conj(L11.'*S1),-conj(L12.'*S1)];
% A_Cof=[L71.'*S1+z*conj(L81.'*S1),conj(L71.'*S1),zero_list,zero_list,...
%     -(L11.'*S1+z*conj(L31.'*S1)),-(L12.'*S1+z*conj(L32.'*S1)),-conj(L11.'*S1),-conj(L12.'*S1)];
B_Cof=[-L2.'*S0-z*conj(L4.'*S0),-conj(L2.'*S0)];
end
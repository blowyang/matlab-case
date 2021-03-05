function B_Cof=B_Matrix_Cof_Three(N,M,R1,R2,x1,x2,a,z)

S0=zeros(M,2*M);
for m=1:M
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i;
end

[~,~,L2,~,~,L4,~,~,~,~,~,~,~]=Get_L_List(N,M,R1,R2,x1,x2,a,z);


B_Cof=[-L2.'*S0-z*conj(L4.'*S0),-conj(L2.'*S0)];
end
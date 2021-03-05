function phi_m_Diff=Get_Phi_Matrix_Diff(a,R,Xc,Xd,z)
N=size(Xc,1)/2;
M=size(Xd,1)/2;
S1=zeros(N,2*N);
S0=zeros(M,2*M);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end
for m=1:M
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i; 
end

L3=zeros(N,1);
L4=zeros(M,1);

for n=1:N
    L3(n,1)=-n/R*(R/z)^(n+1);    
end
for m=1:M
    L4(m,1)=Get_P_Diff(m,a,z);
end   
phi_m_Diff=L3.'*S1*Xc+L4.'*S0*Xd;

end
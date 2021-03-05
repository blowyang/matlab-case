function psi_m=Get_Psi_Matrix(a,R,Xe,Xf,z)
N=size(Xe,1)/2;
M=size(Xf,1)/2;
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

L1=zeros(N,1);
L2=zeros(M,1);

for n=1:N
    L1(n,1)=(R/z)^n;    
end
for m=1:M
    L2(m,1)=Get_P(m,a,z);
end   
psi_m=L1.'*S1*Xe+L2.'*S0*Xf;

end
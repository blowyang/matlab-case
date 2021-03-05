function phi_m_2Diff=Get_Phi_Matrix_2Diff(a,R,Xc,Xd,z)
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

L5=zeros(N,1);
L6=zeros(M,1);

for n=1:N
    L5(n,1)=n*(n+1)/R^2*(R/z)^(n+2);    
end
for m=1:M
    L6(m,1)=Get_P_2Diff(m,a,z);
end   
phi_m_2Diff=L5.'*S1*Xc+L6.'*S0*Xd;

end
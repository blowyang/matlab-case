function S013_Stress_Divided_G0=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,Z)
N=size(Xd,1)/2;
M=size(Xe,1)/2;
S1=zeros(N,2*N);
S2=zeros(M,2*M);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end
for m=1:M
    S2(m,2*m-1)=1;
    S2(m,2*m)=1i; 
end

L1=zeros(N,1);
L2=zeros(M,1);

for n=1:N
    L1(n,1)=n*(R1/Z)^(n+1);
end
for m=1:M
    L2(m,1)=Get_P_Diff(m,a,Z);
end 
S013_Stress_Divided_G0=imag(-1/R1*Xd.'*S1.'*L1+Xe.'*S2.'*L2);
end

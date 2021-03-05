function f0=Get_f0(a,R1,R2,x1,x2,Xcd,Xe,z)
N=size(Xcd,1)/4;
M=size(Xe,1)/2;
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
L2=zeros(N,1);
L3=zeros(M,1);

for n=1:N
    L1(n,1)=(R1/(z-x1))^n;
    L2(n,1)=(R2/(z-x2))^n;
end
for m=1:M
    L3(m,1)=Get_P(m,a,z);
end   
f0=[L1.'*S1,L2.'*S1]*Xcd+L3.'*S0*Xe;
end
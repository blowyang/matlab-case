function [u1,u2]=Get_Dispalcement(a,R,Gm,km,Xce,Xdf,z)
N=size(Xce,1)/4;
M=size(Xdf,1)/4;
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

[L1,L2,L3,L4,~,~]=Get_L_List(N,M,a,R,z);

temp=([km*L1.'*S1-z*conj(L3.'*S1),-conj(L1.'*S1)]*Xce+[km*L2.'*S0-z*conj(L4.'*S0),-conj(L2.'*S0)]*Xdf)/2/Gm;
u1=real(temp);
u2=imag(temp);
end

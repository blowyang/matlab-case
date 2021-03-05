function [Sm22,Sm11,Sm12]=Get_Matrix_Stress(a,R,Xce,Xdf,z)
N=size(Xce,1)/4;
M=size(Xdf,1)/4;

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
[~,~,L3,L4,L5,L6,~,~]=Get_L_List(N,M,R,a,z);

temp=[L3.'*S1+conj(L3.'*S1)+conj(z)*L5.'*S1,conj(L3.'*S1)]*Xce+...
    [L4.'*S0+conj(L4.'*S0)+conj(z)*L6.'*S0,conj(L4.'*S0)]*Xdf;
Sm22=real(temp);
Sm12=imag(temp);
temp=[L3.'*S1+conj(L3.'*S1)-conj(z)*L5.'*S1,-conj(L3.'*S1)]*Xce+...
    [L4.'*S0+conj(L4.'*S0)-conj(z)*L6.'*S0,-conj(L4.'*S0)]*Xdf;
Sm11=real(temp);

end
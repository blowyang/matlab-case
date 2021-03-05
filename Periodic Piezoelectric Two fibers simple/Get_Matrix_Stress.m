function [Sm22,Sm11,Sm12]=Get_Matrix_Stress(R1,R2,x1,x2,a,Xce,Xdf,z)
N=size(Xce,1)/8;

S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i;
end
S0=zeros(N,2*N);
for m=1:N
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i; 
end
[~,~,~,L31,L32,L4,L51,L52,L6,~,~,~,~]=Get_L_List(N,R1,R2,x1,x2,a,z);
temp=[L31.'*S1+conj(L31.'*S1)+conj(z)*L51.'*S1,L32.'*S1+conj(L32.'*S1)+conj(z)*L52.'*S1,conj(L31.'*S1),conj(L32.'*S1)]*Xce+...
    [L4.'*S0+conj(L4.'*S0)+conj(z)*L6.'*S0,conj(L4.'*S0)]*Xdf;
Sm22=real(temp);
Sm12=imag(temp);
temp=[L31.'*S1+conj(L31.'*S1)-conj(z)*L51.'*S1,L32.'*S1+conj(L32.'*S1)-conj(z)*L52.'*S1,-conj(L31.'*S1),-conj(L32.'*S1)]*Xce+...
    [L4.'*S0+conj(L4.'*S0)-conj(z)*L6.'*S0,-conj(L4.'*S0)]*Xdf;
Sm11=real(temp);

end
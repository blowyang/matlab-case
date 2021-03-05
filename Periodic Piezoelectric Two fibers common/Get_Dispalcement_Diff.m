function [u1_Diff,u2_Diff]=Get_Dispalcement_Diff(R1,R2,x1,x2,a,km,Xce,Xdf,z)
N=size(Xce,1)/8;
S1=zeros(N,2*N);
S0=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end
for m=1:N
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i; 
end

[~,~,~,L31,L32,L4,L51,L52,L6,~,~,~,~]=Get_L_List(N,R1,R2,x1,x2,a,z);

temp=[km*L31.'*S1-conj(L31.'*S1)-z*conj(L51.'*S1),km*L32.'*S1-conj(L32.'*S1)-z*conj(L52.'*S1),...
    -conj(L31.'*S1),-conj(L32.'*S1)]*Xce+...
    [km*L4.'*S0-conj(L4.'*S0)-z*conj(L6.'*S0),-conj(L4.'*S0)]*Xdf;
u1_Diff=real(temp);
u2_Diff=imag(temp);
end

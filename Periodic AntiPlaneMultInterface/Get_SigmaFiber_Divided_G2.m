function [Sigma_zs_Divided_G2,Sigma_zr_Divided_G2,Sigma_23_Divided_G2,Sigma_13_Divided_G2]=Get_SigmaFiber_Divided_G2(Xa,R2,Z)
%获得fiber上的环向和轴向应力
N=size(Xa,1)/2;
S1=zeros(N,2*N);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end

L1=zeros(N,1);

for n=1:N
    L1(n,1)=n*(Z/R2)^(n-1);
end

Sigma_23_Divided_G2=real(1/R2*Xa.'*S1.'*L1);
Sigma_13_Divided_G2=imag(1/R2*Xa.'*S1.'*L1);

Sigma_zs_Divided_G2=real(exp(1i*angle(Z))*(Sigma_23_Divided_G2+Sigma_13_Divided_G2*1i));
Sigma_zr_Divided_G2=imag(exp(1i*angle(Z))*(Sigma_23_Divided_G2+Sigma_13_Divided_G2*1i));

end
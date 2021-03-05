function res=dpsi(S11,S22,S12,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)
% B1=(S11+S22)/4;
% C1=0;
B2=(S22-S11)/2;
C2=S12;

res=(B2+1i*C2)+conj(gama)/(Z-Zd)+gama*conj(Zd)/(Z-Zd)^2;
for k=1:1:N
    ck=X(3*N+k);
    dk=X(4*N+k);
    jk=X(5*N+k);
    res=res-k*ck/R1*(R1/(Z-Z1))^(k+1)-k*dk/R2*(R2/(Z-Z2))^(k+1)-k*jk/R3*(R3/(Z-Z3))^(k+1);
end
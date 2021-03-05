function res=dfi(S11,S22,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)
B1=(S11+S22)/4;
C1=0;
% B2=(S22-S11)/2;
% C2=S12;
res=(B1+1i*C1)+gama/(Z-Zd);
for k=1:1:N
    ak=X(k);
    bk=X(N+k);
    ik=X(2*N+k);
    res=res-k*ak/R1*(R1/(Z-Z1))^(k+1)-k*bk/R2*(R2/(Z-Z2))^(k+1)-k*ik/R3*(R3/(Z-Z3))^(k+1);
end
function res=ddfi(N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)
res=-gama/(Z-Zd)^2;
for k=1:1:N
    ak=X(k);
    bk=X(N+k);
    ik=X(2*N+k);
    res=res+k*(k+1)*ak/R1/R1*(R1/(Z-Z1))^(k+2)+k*(k+1)*bk/R2/R2*(R2/(Z-Z2))^(k+2)+k*(k+1)*ik/R3/R3*(R3/(Z-Z3))^(k+2);
end
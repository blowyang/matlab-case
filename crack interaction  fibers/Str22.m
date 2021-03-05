function res=Str22(S11,S22,S12,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)

res=real(dfi(S11,S22,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)+conj(dfi(S11,S22,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama))+conj(Z)*ddfi(N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama)+dpsi(S11,S22,S12,N,R1,R2,R3,Z1,Z2,Z3,X,Z,Zd,gama));

function res=FX(N,R1,R2,Z1,Z2,X,Zd,bx,by)

S22=real(dfi(N,R1,R2,Z1,Z2,X,Zd)+conj(dfi(N,R1,R2,Z1,Z2,X,Zd))+conj(Zd)*ddfi(N,R1,R2,Z1,Z2,X,Zd)+dpsi(N,R1,R2,Z1,Z2,X,Zd));
%S11=real(dfi(N,R1,R2,Z1,Z2,X,Zd)+conj(dfi(N,R1,R2,Z1,Z2,X,Zd))-conj(Zd)*ddfi(N,R1,R2,Z1,Z2,X,Zd)-dpsi(N,R1,R2,Z1,Z2,X,Zd));
S12=imag(dfi(N,R1,R2,Z1,Z2,X,Zd)+conj(dfi(N,R1,R2,Z1,Z2,X,Zd))+conj(Zd)*ddfi(N,R1,R2,Z1,Z2,X,Zd)+dpsi(N,R1,R2,Z1,Z2,X,Zd));
res=S12*bx+S22*by;

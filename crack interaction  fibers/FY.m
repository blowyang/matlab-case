function res=FY(N,R1,R2,Z1,Z2,RX,IX,Zd,bx,by)

%S22=real(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd)+conj(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd))+conj(Zd)*ddfi(N,R1,R2,Z1,Z2,RX,IX,Zd)+dpsi(N,R1,R2,Z1,Z2,RX,IX,Zd));
S11=real(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd)+conj(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd))-conj(Zd)*ddfi(N,R1,R2,Z1,Z2,RX,IX,Zd)-dpsi(N,R1,R2,Z1,Z2,RX,IX,Zd));
S12=imag(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd)+conj(dfi(N,R1,R2,Z1,Z2,RX,IX,Zd))+conj(Zd)*ddfi(N,R1,R2,Z1,Z2,RX,IX,Zd)+dpsi(N,R1,R2,Z1,Z2,RX,IX,Zd));
res=-S11*bx-S12*by;
